(ns
  signature-invariants.signature
  "Methods to algebraically work with the signature."
  (:require [hopf-algebra.shuffle-algebra :as sa]
            [hopf-algebra.linear-combination :as lc]
            [hopf-algebra.linear-combination-plus :as lcp]
            [hopf-algebra.hopf-algebra :refer [product]]; coproduct antipode to-str to-latex]]
            [clojure.math.numeric-tower :refer [expt]]
            [spyscope.core]
            ))

(defn zeroes [n]
  (vec (take n (repeat 0))))

(defn fact
  "n !"
  [n]
  (apply * (range 1 (inc n))))

(defn n-choose-k [n k]
  (/ (fact n) (* (fact k) (fact (- n k)))))

(defn- as-vector-get-position [m dim order]
  (loop [i 0,
         result 0]
    (if (>= i (count m))
      result
      (recur
        (inc i)
        (+ result (* (expt dim (- order i 1)) (nth m i)))))))

(defn as-vector
  "Given a linear combination `lc` of ShuffleWord/ConcatWord
   in dimension `dim` of order `order`,
   return the corresponding vector in \\R^{dim * order}."
  [lc dim order]
  (loop [result (zeroes (expt dim order)),
         remaining lc]
    (if (empty? remaining)
      result
      (recur
        (assoc result
               (as-vector-get-position (:content (key (first remaining))) dim order)
               (val (first remaining)))
        (rest remaining)))))

(defn x-dy
  "< S(X), x-dy > := <S(X), x> d<S(X), y>"
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
  [lc-1 lc-2]
  (lc/apply-bilinear-function
    (fn [x y] (lc/+ (x-dy x y) (lc/* -1 (x-dy y x))))
    lc-1 lc-2))


(defn- subword [coll n m]
  (take
    (- m n)
    (drop n coll)))

(defn- all-splits
  "Splits list `ell` into all possible nonempty lists.

  [:a :b] => [ [:a :b] ] [[:a] [:b] ]"
  [ell]
  (for [ss (clojure.math.combinatorics/subsets (range 0 (dec (count ell))))]

    (loop [current 0
           remaining-sss (sort ss)
           result []
           ]
      (if (empty? remaining-sss)
        (conj result
              (subword ell current (count ell)))
        (recur (inc (first remaining-sss))
               (rest remaining-sss)
                (conj result
                  (subword ell current (inc (first remaining-sss)))))))))

(defn dual-logarithm
 " < dual-logarithm(psi), x > = < psi, log(x) >,
  compare Reutenauer - Free Lie algebras (3.2.3)"
  [sw]
  (loop [splits (all-splits (:content sw))]
    (apply lc/+
      (map
        (fn [s]
          (let [n (count s)]
            (lc/*
              (/ (expt -1 (inc n)) n)
              (apply lc/*
                     (map (fn [w] {(sa/->ShuffleWord w) 1}) s)))))
        splits))))
  

(defn lc-dual-logarithm [lc] ; XXX name
  (lc/apply-linear-function dual-logarithm lc))

;; OLD IMPLEMENTATION
;;(defn- de-concat
;;  "Split shuffle-word up into m, possibly empty, words."
;;  [m shuffle-word]
;;  (if (= 2 m)
;;    (sa/deconcatenation (:content shuffle-word))
;;    (lc/apply-linear-function
;;      (fn [ [x & r] ]
;;        (lc/otimes
;;          (sa/deconcatenation (:content x))
;;          (lc/lift r)))
;;    (de-concat (dec m) shuffle-word))))
;;
;;(defn- shuffle-de-concat 
;;  "= \\sum_{p_1..p_m = shuffle-word} p_1 \\shuffle .. \\shuffle p_m "
;;  [m shuffle-word]
;;  (if (= 1 m)
;;    (lc/lift shuffle-word)
;;    (lc/apply-linear-function
;;      (fn [tensor]
;;        (apply
;;          lc/*
;;          (map lc/lift tensor)))
;;      (de-concat m shuffle-word))))
;;
;;(defn- c [N m]
;;  (reduce +
;;          (map
;;            (fn [r]
;;              (*
;;                (/ (expt -1 (+ (* 2 r) (- 1 m))) r)
;;                (n-choose-k r m)))
;;            (range m (inc N)))))
;;
;;(defn dual-logarithm
;;  " < dual-logarithm(psi), x > = < psi, log(x) > "
;;  ([shuffle-word] (dual-logarithm shuffle-word (count (:content shuffle-word))))
;;  ([shuffle-word NN]
;;  (let [N NN]
;;    (reduce
;;      lc/+
;;      (map
;;        (fn [m]
;;          (lc/*
;;            (c N m)
;;            (shuffle-de-concat m shuffle-word)))
;;        (range 1 (inc N)))))))




(defn project
  "Projection to \\le N"
  [x N]
  (lc/apply-linear-function
    (fn [x] 
      (if (> (count (:content x)) N)
        {}
        (lc/lift x)))
    x))

(defn- exponential-helper [lie n N]
  (if (= 0 n)
    {(sa/->ConcatWord []) 1}
    (if (= 1 n)
      lie
      (reduce
        (fn [a b]
          (project (lc/* a b) N))
        (repeat n lie)))))

; XXX this does not allow for symbolic coefficients ..
(defn exponential
  "Exponential of `lie` in the concatenation algebra, up to level \\le `N`."
  [lie N]
  (reduce
    lc/+
    (map
      (fn [n]
        (lc/* (/ 1 (fact n))
                        (exponential-helper lie n N)))
      (range 0 (inc N)))))



(defn- project-symbolic [x N]
  (lcp/apply-linear-function
    (fn [x] ; projection to \le N
      (if (> (count (:content x)) N)
        {}
        {x 1}
        ;{x {(sa/->ConcatWord []) 1}}
        ))
    x))
(defn- exponential-symbolic-helper [lie n N unit]
  (if (= 0 n)
    {(sa/->ConcatWord []) unit}
     ;{(sa/->ConcatWord []) 1}} ; XXX
    (if (= 1 n)
      lie
      (reduce
        (fn [a b]
          (project-symbolic (lcp/* a b) N))
        (repeat n lie)))))
(defn exponential-symbolic
  "lie: linear combination plus of lie elements
   N: project to level <= N
   unit: the unit in the ring of coefficients"
  ([lie unit]
   (let [N (last (sort (map (fn [ [k v] ] (count (:content k))) lie)))]
     (exponential-symbolic lie unit N)))
  ([lie unit N]
  (reduce
    lcp/+
    (map
      (fn [n] (lcp/* (/ 1 (fact n))
                     ;#spy/d
                     (exponential-symbolic-helper lie n N unit)))
      (range 0 (inc N))))))


;(defn inverse-symbolic [g N]
;  (assert (= {(sa/->ConcatWord []) 1} (get g (sa/->ConcatWord []))))
;  (let [z (lc/lc-plus-add g
;                          {(sa/->ConcatWord []) {(sa/->ConcatWord []) -1}})]
;  (reduce
;    lc/lc-plus-add
;    (map
;      (fn [n] (lc/lc-plus-multiply-scalar (exp -1 n)
;                                          ;#spy/d
;                                          (exponential-symbolic-helper z n N)))
;      (range 0 (inc N))))))
