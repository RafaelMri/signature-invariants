(ns
  ^{:doc
        ".."}
  signature-invariants.signature
  (:require [hopf-algebra.shuffle-algebra :as sa]
            [hopf-algebra.linear-combination :as lc]
            [hopf-algebra.hopf-algebra :refer [product]]; coproduct antipode to-str to-latex]]
            [spyscope.core]
            ))

(defn zeroes [n]
  (vec (take n (repeat 0))))

(defn exp
  "x^n for n a natural number."
  [x n]
  (reduce * (repeat n x)))

(defn fact
  "n ! for n a natural number."
  [n] (apply * (range 1 (inc n))))

(defn n-choose-k [n k]
  (/ (fact n) (* (fact k) (fact (- n k)))))

(defn- as-vector-get-position [m dim order]
  (loop [i 0,
         result 0]
    (if (>= i (count m))
      result
      (recur
        (inc i)
        (+ result (* (exp dim (- order i 1)) (nth m i)))))))

(defn as-vector
  "Given a linear combination `lc` of ShuffleWord/ConcatWord
   in dimension `dim` of order `order`,
   return the corresponding vector in \\R^{dim * order}."
  [lc dim order]
  (loop [result (zeroes (exp dim order)),
         remaining lc]
    (if (empty? remaining)
      result
      (recur
        (assoc result
               (as-vector-get-position (:content (key (first remaining))) dim order)
               (val (first remaining)))
        (rest remaining)))))

(defn x-dy
  "<S(X), x> d<S(X), y> = < S(X), x-dy >."
  [x-w y-w]
  (let [x (:content x-w)
        y (:content y-w)]
  (if (or (empty? x) (empty? y))
    {}
    (let [b-beginning (take (dec (count y)) y),
          integrands (product (sa/->ShuffleWord x) (sa/->ShuffleWord b-beginning))]
      (into
        {}
        (map
          (fn [ [k v] ] [ (sa/->ShuffleWord (vec (concat (:content k) [(last y)]))) v ])
          integrands))))))

(defn area-between
  [a b]
  (lc/lc-apply-bilinear-function
    (fn [x y] (lc/lc-add (x-dy x y) (lc/lc-multiply -1 (x-dy y x))))
    a b))


(defn- de-concat
  [m shuffle-word]
  (if (= 2 m)
    (sa/deconcatenation (:content shuffle-word))
    (lc/lc-apply-linear-function
      (fn [ [x & r] ]
        (lc/lc-otimes
          (sa/deconcatenation (:content x))
          (lc/lc-lift r)))
    (de-concat (dec m) shuffle-word)
      )))

(defn- de-shuffle
  "= \\sum_{p_1..p_m} p_1 \\shuffle .. \\shuffle p_m "
  [m shuffle-word]
  (if (= 1 m)
    (lc/lc-lift shuffle-word)
    (lc/lc-apply-linear-function
      (fn [tensor]
        (apply
          lc/lc-multiply
          (map lc/lc-lift tensor)))
      (de-concat m shuffle-word))))

(defn- c [N m]
  (reduce +
          (map
            (fn [r]
              (*
                (/ (exp -1 (+ (* 2 r) (- 1 m))) r)
                (n-choose-k r m)))
            (range m (inc N)))))

(defn dual-exponential
  " < dual-exponential(psi), x > = < psi, log(x) > "
  ([shuffle-word] (dual-exponential shuffle-word (count (:content shuffle-word))))
  ([shuffle-word NN]
  (let [N NN]
    (reduce
      lc/lc-add
      (map
        (fn [m]
          (lc/lc-multiply
            (c N m)
            (de-shuffle m shuffle-word)))
        (range 1 (inc N)))))))

(defn lc-dual-exponential [lc]
  (lc/lc-apply-linear-function
    dual-exponential
    lc))


(defn project
  "Projection to \\le N"
  [x N]
  (lc/lc-apply-linear-function
    (fn [x] 
      (if (> (count (:content x)) N)
        {}
        (lc/lc-lift x)))
    x))

(defn- exponential-helper [lie n N]
  (if (= 0 n)
    {(sa/->ConcatWord []) 1}
    (if (= 1 n)
      lie
      (reduce
        (fn [a b]
          (project (lc/lc-multiply a b) N))
        (repeat n lie)))))

(defn exponential
  "Exponential of `lie` in the concatenation algebra, up to level \\le `N`."
  [lie N]
  (reduce
    lc/lc-add
    (map
      (fn [n]
        (lc/lc-multiply (/ 1 (fact n))
                        (exponential-helper lie n N)))
      (range 0 (inc N)))))
