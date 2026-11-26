/-!
# The table of the base-width Gaussian sampler

`rust_implem/src/sampling.rs` samples the centered discrete Gaussian of width
sigma = 3.2 from a table. It draws a uniform 63-bit integer `u`, takes the
number of thresholds `CDT[k] <= u` as the magnitude and gives it a uniform sign.
Magnitude k therefore has probability `(CDT[k] - CDT[k-1]) / 2^63`, with
`CDT[-1] = 0` and `CDT[29] = 2^63`.

This file

  1. computes the ideal distribution of the magnitude in fixed point,
  2. builds the table from it: every magnitude k >= 1 has its ideal probability
     rounded to the nearest multiple of 2^-63, and zero takes what remains,
  3. reads the Rust source and checks that its constant is this table,
  4. checks that the table gives no tail more probability than the ideal
     distribution does, and that zero has to take what remains for that,
  5. measures how far the sampled distribution is from the ideal one, and does
     the same for other ways of filling a 63-bit table,
  6. bounds what any table of this precision can reach.

The measures compare the sampled distribution Q of the magnitude with the ideal
one P (the sign is uniform and independent in both, so the signed distributions
give the same values):

  * statistical distance: half the sum of `|Q(k) - P(k)|`,
  * Renyi divergence of order 2: the sum of `Q(k)^2 / P(k)`,
  * Renyi divergence of order infinity: the largest `Q(k) / P(k)`.

## What the theorems are

Every theorem below is a closed statement that Lean decides by evaluation:
`native_decide` compiles the statement and runs it, so it relies on the compiler
of Lean as well as on its kernel.

The weights `exp(-k^2 / (2 sigma^2))` are computed with 512 fractional bits by a
Taylor series whose steps round down. They are approximations: nothing here
relates them to the exponential of the real numbers. `reference_is_stable`
shows that computing them with 768 bits changes the reference by less than
2^-508, far below anything that the statements depend on. The unit test
`cdt_matches_exact_distribution` of the Rust crate computes the table once more,
in another way.

## Running

    cd tools/lean && lean ./cdt_table.lean

prints `cdt_table_out.txt`. The version of Lean is in `lean-toolchain`. Lean
needs the `./` to find the Rust source.
-/

/-! ## The ideal distribution -/

/-- sigma^2 = 1024 / 100. -/
def SIGMA2_NUM : Nat := 1024
def SIGMA2_DEN : Nat := 100

/-- Fractional bits of the fixed-point numbers. -/
def PREC : Nat := 512

/-- `exp(-x)` in fixed point: for `X = x * 2^P`, close to `exp(-x) * 2^P`. -/
def expNeg (X : Nat) (P : Nat := PREC) : Nat := Id.run do
  let one : Nat := 1 <<< P
  -- x < 2^m, so y = x / 2^m is below 1.
  let m : Nat := Nat.log2 (X >>> P) + 1
  let y : Nat := X >>> m
  -- exp(y) by its Taylor series.
  let mut term : Nat := one
  let mut sum : Nat := one
  for i in [1:200] do
    term := ((term * y) >>> P) / i
    sum := sum + term
  -- exp(-x) = (1 / exp(y)) ^ (2^m).
  let mut v : Nat := (one * one) / sum
  for _ in [0:m] do
    v := (v * v) >>> P
  return v

/-- Cumulative distribution of the magnitude over the magnitudes below `n`:
    entry k is `P(|x| <= k) * 2^P`, and the last entry is exactly `2^P`. -/
def idealCumulative (n : Nat) (P : Nat := PREC) : Array Nat := Id.run do
  -- Weight of a magnitude: rho(0), then 2 rho(k), with rho(k) = exp(-k^2 / (2 sigma^2)).
  let mut w : Array Nat := #[]
  for k in [0:n] do
    let rho := expNeg (((k * k * SIGMA2_DEN) <<< P) / (2 * SIGMA2_NUM)) P
    w := w.push (if k == 0 then rho else 2 * rho)
  let norm := w.foldl (fun a b => a + b) 0
  let mut cumulative : Array Nat := #[]
  let mut acc : Nat := 0
  for x in w do
    acc := acc + x
    cumulative := cumulative.push ((acc <<< P) / norm)
  return cumulative

/-- The reference: magnitudes 0 to 49. Those left out have a probability below 2^-178. -/
def reference : Array Nat := idealCumulative 50

/-- Consecutive differences: the probabilities behind cumulative values. -/
def differences (cumulative : Array Nat) : Array Nat := Id.run do
  let mut d : Array Nat := #[]
  let mut prev : Nat := 0
  for x in cumulative do
    d := d.push (x - prev)
    prev := x
  return d

/-- Ideal probability of each magnitude, in units of `2^-PREC`. -/
def ideal : Array Nat := differences reference

/-- No measure below divides by zero. -/
theorem ideal_is_positive : ideal.all (fun p => p > 0) = true := by native_decide

/-- With 768 fractional bits the reference is the same to within 2^-508. -/
theorem reference_is_stable :
    (List.range 50).all (fun k =>
      let fine := (idealCumulative 50 768)[k]! >>> 256
      let coarse := reference[k]!
      (if fine >= coarse then fine - coarse else coarse - fine) < 2 ^ 4) = true := by
  native_decide

/-! ## Tables

A table is the array of its thresholds, in units of `2^-63`. -/

/-- `2^-63` in units of `2^-PREC`. -/
def UNIT : Nat := 1 <<< (PREC - 63)

/-- Nearest multiple of `2^-63`. -/
def nearest (x : Nat) : Nat := (x + UNIT / 2) / UNIT

/-- Largest multiple of `2^-63` that is not above `x`. -/
def below (x : Nat) : Nat := x / UNIT

/-- Thresholds for the probabilities `q` of the magnitudes 0, 1, ... (in units of
    `2^-63`, `2^63` in total): the cumulative sums that are below `2^63`. -/
def thresholds (q : Array Nat) : Array Nat := Id.run do
  let mut t : Array Nat := #[]
  let mut acc : Nat := 0
  for x in q do
    acc := acc + x
    if acc < 2 ^ 63 then t := t.push acc
  return t

/-- The construction of the table of the sampler: every magnitude k >= 1 has its
    ideal probability rounded to the nearest multiple of `2^-63`, and zero takes
    what remains. -/
def roundedProbabilities : Array Nat :=
  let q := (ideal.extract 1 ideal.size).map nearest
  thresholds (#[2 ^ 63 - q.foldl (fun a b => a + b) 0] ++ q)

/-- Every probability rounded to the nearest multiple of `2^-63`, that of zero
    included. -/
def rounded : Array Nat := ideal.map nearest

/-- They add up to `2^63 - 1`: one unit is left. -/
theorem rounding_leaves_one_unit : rounded.foldl (fun a b => a + b) 0 + 1 = 2 ^ 63 := by
  native_decide

/-- The table in which magnitude `k` takes that unit. -/
def unitOn (k : Nat) : Array Nat := thresholds (rounded.modify k (fun q => q + 1))

/-- In the table of the sampler, zero takes it. -/
theorem zero_takes_the_unit : roundedProbabilities = unitOn 0 := by native_decide

/-- Cumulative values rounded to nearest. The one of magnitude 29 rounds to
    `2^63`, which no draw reaches. -/
def roundedCumulative : Array Nat := (reference.map nearest).filter (fun t => t < 2 ^ 63)

/-- Cumulative values rounded down, up to magnitude 29: beyond it they are all
    `2^63 - 1`. -/
def flooredCumulative : Array Nat := (reference.extract 0 30).map below

/-- Cumulative values of the distribution cut at magnitude 29, rounded down. The
    last one is `2^63`, which no draw reaches. -/
def flooredCut : Array Nat := ((idealCumulative 30).map below).filter (fun t => t < 2 ^ 63)

/-- Weights in double precision for the magnitudes 0 to 39, divided by their sum,
    accumulated, scaled by `2^63` and rounded down (with the `exp` of glibc:
    another one can change the last bits). Just below 1 the doubles are `2^-53`
    apart and an increment below `2^-54` is lost: the entries of the magnitudes
    27 to 39 are equal. -/
def doublePrecision : Array Nat := #[
  1149872835429265536, 3340023666152831488, 5231742854224523264, 6713673034491316224,
  7766573326200194048, 8445050402542554112, 8841576285654609920, 9051758678878183424,
  9152802451769413632, 9196859074767744000, 9214281206174001152, 9220529764022705152,
  9222562339745869824, 9223161995634593792, 9223322447917708288, 9223361386320108544,
  9223369956674608128, 9223371667508610048, 9223371977254383616, 9223372028116137984,
  9223372035690843136, 9223372036713967616, 9223372036839305216, 9223372036853230592,
  9223372036854633472, 9223372036854761472, 9223372036854771712, 9223372036854772736,
  9223372036854772736, 9223372036854772736, 9223372036854772736, 9223372036854772736,
  9223372036854772736, 9223372036854772736, 9223372036854772736, 9223372036854772736,
  9223372036854772736, 9223372036854772736, 9223372036854772736, 9223372036854772736]

/-! ## The table of the Rust source -/

def rustSource : String := include_str "../../rust_implem/src/sampling.rs"

/-- The entries of `const CDT` in the Rust source. -/
def table : Array Nat :=
  let constant := (rustSource.splitOn "\nconst CDT: [u64; ")[1]!
  let entries := ((constant.splitOn "= [")[1]!.splitOn "];")[0]!
  let digits := String.ofList (entries.toList.filter (fun c => c.isDigit || c == ','))
  ((digits.splitOn ",").filter (fun s => s != "")).toArray.map (fun s => s.toNat!)

/-- The sampler uses the table built above. -/
theorem table_is_rounded_probabilities : table = roundedProbabilities := by native_decide

/-! ## Tails -/

/-- Whether no threshold of `t` is below the ideal cumulative value: `t` then
    gives no tail `|x| > k` more probability than the discrete Gaussian does.
    (Cut at magnitude 49, the reference is if anything above the cumulative
    values of the whole distribution.) -/
def noHeavierTail (t : Array Nat) : Bool :=
  (List.range t.size).all (fun k => t[k]! * UNIT >= reference[k]!)

theorem no_tail_is_heavier : noHeavierTail table = true := by native_decide

/-- Only zero can take the unit that rounding leaves: on any other magnitude,
    some threshold is below the ideal cumulative value. -/
theorem only_zero_keeps_the_tails :
    (List.range 50).all (fun k => k == 0 || !noHeavierTail (unitOn k)) = true := by
  native_decide

/-! ## Measures -/

/-- Probabilities that a table gives to the magnitudes 0, 1, ..., in units of
    `2^-PREC`: `t[k] - t[k-1]`, and for the largest magnitude what remains up to
    `2^63`. -/
def sampled (t : Array Nat) : Array Nat :=
  (differences (t.push (2 ^ 63))).map (fun q => q * UNIT)

/-- A fraction: numerator and denominator. -/
abbrev Ratio := Prod Nat Nat

/-- Statistical distance between the distribution that `t` samples and the ideal one. -/
def statisticalDistance (t : Array Nat) : Ratio := Id.run do
  let q := sampled t
  let mut sum : Nat := 0
  for k in [0:max q.size ideal.size] do
    let a := q.getD k 0
    let b := ideal.getD k 0
    sum := sum + (if a >= b then a - b else b - a)
  return (sum, 2 * 2 ^ PREC)

/-- Renyi divergence of order 2 of the sampled distribution from the ideal one, minus 1. -/
def renyi2Excess (t : Array Nat) : Ratio := Id.run do
  let q := sampled t
  let mut sum : Nat := 0
  for k in [0:q.size] do
    sum := sum + q[k]! * q[k]! / ideal[k]!
  return (sum - 2 ^ PREC, 2 ^ PREC)

/-- Renyi divergence of order infinity of the sampled distribution from the ideal
    one, minus 1: the largest ratio of a sampled probability to the ideal one. -/
def renyiInfExcess (t : Array Nat) : Ratio := Id.run do
  let q := sampled t
  let mut largest : Ratio := (1, 1)
  for k in [0:q.size] do
    if q[k]! * largest.2 > largest.1 * ideal[k]! then largest := (q[k]!, ideal[k]!)
  return (largest.1 - largest.2, largest.2)

/-- `r < s`. -/
def less (r s : Ratio) : Bool := r.1 * s.2 < s.1 * r.2

/-- `r < 2^(-a/100)`. -/
def belowTwoToMinus (r : Ratio) (a : Nat) : Bool := r.1 ^ 100 * 2 ^ a < r.2 ^ 100

/-- 2^-60.97 <= statistical distance < 2^-60.96. -/
theorem statistical_distance :
    belowTwoToMinus (statisticalDistance table) 6096 = true /\
    belowTwoToMinus (statisticalDistance table) 6097 = false := by native_decide

/-- 2.64e-20 < Renyi divergence of order 2, minus 1 < 2.65e-20. -/
theorem renyi_divergence_of_order_2 :
    less (264, 10 ^ 22) (renyi2Excess table) = true /\
    less (renyi2Excess table) (265, 10 ^ 22) = true := by native_decide

/-- 1.18e-5 < Renyi divergence of order infinity, minus 1 < 1.19e-5. -/
theorem renyi_divergence_of_order_infinity :
    less (118, 10 ^ 7) (renyiInfExcess table) = true /\
    less (renyiInfExcess table) (119, 10 ^ 7) = true := by native_decide

/-- On each measure the table of the sampler is closer to the ideal distribution
    than the tables that are built from cumulative values. -/
theorem closer_than_tables_of_cumulative_values :
    [roundedCumulative, flooredCumulative, flooredCut, doublePrecision].all (fun other =>
      less (statisticalDistance table) (statisticalDistance other) &&
      less (renyi2Excess table) (renyi2Excess other) &&
      less (renyiInfExcess table) (renyiInfExcess other)) = true := by native_decide

/-- Where the unit that rounding leaves would bring the table closest: on 5 for
    the Renyi divergence of order 2, on 28 for the statistical distance. -/
theorem closest_places_for_the_unit :
    (List.range 50).all (fun k =>
      (k == 5 || less (renyi2Excess (unitOn 5)) (renyi2Excess (unitOn k))) &&
      (k == 28 || less (statisticalDistance (unitOn 28)) (statisticalDistance (unitOn k))))
      = true := by native_decide

/-! ## What 63 bits allow

Whatever its thresholds, a table gives every magnitude a probability Q(k) that
is a multiple of `2^-63`, at least as far from the ideal probability P(k) as the
nearest multiple is: half the sum of these gaps is a lower bound of the
statistical distance. Both distributions sum to one, so the Renyi divergence of
order 2 minus 1 is the sum of `(Q(k) - P(k))^2 / P(k)`, which is at least the
sum of `gap^2 / P(k)`. -/

/-- Distance of each ideal probability from the nearest multiple of `2^-63`. -/
def gaps : Array Nat := ideal.map (fun p => min (p % UNIT) (UNIT - p % UNIT))

/-- No table has a smaller statistical distance. -/
def smallestDistance : Ratio := (gaps.foldl (fun a b => a + b) 0, 2 * 2 ^ PREC)

/-- No table has a smaller Renyi divergence of order 2, minus 1. -/
def smallestRenyi2Excess : Ratio := Id.run do
  let mut sum : Nat := 0
  for k in [0:ideal.size] do
    sum := sum + gaps[k]! * gaps[k]! / ideal[k]!
  return (sum, 2 ^ PREC)

/-- The statistical distance of the table is within 9% of the lower bound. -/
theorem statistical_distance_is_close_to_the_smallest :
    less (statisticalDistance table)
      (109 * smallestDistance.1, 100 * smallestDistance.2) = true := by native_decide

/-- Its Renyi divergence of order 2, minus 1, is within 15 digits of the lower bound. -/
theorem renyi_divergence_of_order_2_is_close_to_the_smallest :
    less (renyi2Excess table)
      ((10 ^ 15 + 1) * smallestRenyi2Excess.1, 10 ^ 15 * smallestRenyi2Excess.2) = true := by
  native_decide

/-! ## Report -/

/-- Three decimal digits. -/
def digits3 (n : Nat) : String := s!"{n / 100}{n / 10 % 10}{n % 10}"

/-- Decimal digits in groups of three, as in the Rust source. -/
partial def grouped (n : Nat) : String :=
  if n < 1000 then toString n else grouped (n / 1000) ++ "_" ++ digits3 (n % 1000)

/-- `n / 100` with two decimals. -/
def hundredth (n : Nat) : String := s!"{n / 100}.{n / 10 % 10}{n % 10}"

/-- `2^(-(a+1)/100) < r <= 2^(-a/100)`. -/
def hundredths (r : Ratio) : Nat := Nat.log2 (r.2 ^ 100 / r.1 ^ 100)

/-- `r` with four significant digits, rounded down. -/
def scientific (r : Ratio) : String := Id.run do
  -- 10^(e - 400) <= r < 10^(e - 399)
  let mut e : Nat := 0
  for i in [0:800] do
    let reached := if i >= 400 then r.1 >= r.2 * 10 ^ (i - 400) else r.1 * 10 ^ (400 - i) >= r.2
    if reached then e := i
  let m :=
    if e >= 400 then r.1 * 1000 / (r.2 * 10 ^ (e - 400)) else r.1 * 1000 * 10 ^ (400 - e) / r.2
  let exponent := if e >= 400 then toString (e - 400) else s!"-{400 - e}"
  return s!"{m / 1000}.{digits3 (m % 1000)}e{exponent}"

def column (s : String) (width : Nat) : String :=
  s ++ String.ofList (List.replicate (width - s.length) ' ')

def row (name : String) (t : Array Nat) : String :=
  let q := sampled t
  let largest := (List.range q.size).foldl (fun m k => if q[k]! > 0 then k else m) 0
  column name 47 ++ column (toString t.size) 9 ++ column (toString largest) 11 ++
    column ("2^-" ++ hundredth (hundredths (statisticalDistance t))) 13 ++
    column (scientific (renyi2Excess t)) 12 ++ scientific (renyiInfExcess t)

#eval do
  IO.println "Table of the base-width Gaussian sampler (const CDT in rust_implem/src/sampling.rs)"
  IO.println ""
  for i in [0:(table.size + 2) / 3] do
    let line := (table.extract (3 * i) (3 * i + 3)).map (fun t => grouped t ++ ",")
    IO.println ("    " ++ String.intercalate " " line.toList)
  IO.println ""
  -- By how much the thresholds exceed the ideal cumulative values, in hundredths of a unit,
  -- rounded down: the largest is printed with one more, as an upper bound.
  let above := (List.range table.size).map
    (fun k => (table[k]! * UNIT - reference[k]!) * 100 / UNIT)
  IO.println ("Every threshold is above the ideal cumulative value, by " ++
    hundredth (above.foldl min above[0]!) ++ " to " ++ hundredth (above.foldl max 0 + 1) ++
    " units of 2^-63:")
  IO.println "no tail |x| > k is heavier than in the discrete Gaussian."
  IO.println ""
  IO.println "Distance of the sampled distribution from the discrete Gaussian, per sample"
  IO.println ""
  IO.println (column "table" 47 ++ column "entries" 9 ++ column "largest" 11 ++
    column "statistical" 13 ++ column "Renyi 2" 12 ++ "Renyi infinity")
  IO.println (column "" 56 ++ column "magnitude" 11 ++ column "distance" 13 ++
    column "minus 1" 12 ++ "minus 1")
  IO.println (row "sampler: probabilities rounded, rest to zero" table)
  IO.println (row "same, the last unit on magnitude 5" (unitOn 5))
  IO.println (row "same, the last unit on magnitude 28" (unitOn 28))
  IO.println (row "cumulative values rounded to nearest" roundedCumulative)
  IO.println (row "cumulative values rounded down" flooredCumulative)
  IO.println (row "cut at magnitude 29, cumulative rounded down" flooredCut)
  IO.println (row "accumulated in double precision" doublePrecision)
  IO.println ""
  IO.println "In the table, the statistical distances are upper bounds and the divergences"
  IO.println "are rounded down."
  IO.println ""
  IO.println "Rounded to nearest, the probabilities add up to 2^63 - 1. Zero takes the unit"
  IO.println "that is left: on any other magnitude, a threshold falls below the ideal"
  IO.println "cumulative value. On 5 the unit gives the smallest Renyi divergence of order 2,"
  IO.println "on 28 the smallest statistical distance."
  IO.println ""
  IO.println ("No table with 63-bit thresholds has a statistical distance below 2^-" ++
    hundredth (hundredths smallestDistance + 1) ++ ",")
  IO.println ("nor a Renyi divergence of order 2 below 1 + " ++
    scientific smallestRenyi2Excess ++ ".")
