(ns
  signature-invariants.lyndon-words
  "Methods to calculate Lyndon words and the Lyndon basis of the free Lie algebra."
  (:require [signature-invariants.young-tableau :as yt]
            [hopf-algebra.linear-combination :as lc]
            [hopf-algebra.shuffle-algebra :as sa]
            [hopf-algebra.hopf-algebra :refer [product]]))

(defn lazy-concat-repeat [coll]
  (lazy-cat coll (lazy-concat-repeat coll)))


(defn remove-trailing [x coll]
    (loop [upto (dec (count coll))]
      (if (or
            (= -1 upto)
            (not (= x (nth coll upto))))
        (take (inc upto) coll)
        (recur (dec upto)))))

(defn lyndon-words
  "Generate all Lyndon words in the letters 0,..,`dim`-1 of length <= length, using Duval's algorithm."
  [dim length]

  (loop [current [-1]
         result []]
    (if (empty? current)
      result
      (let [current (update current (dec (count current)) inc)]
        (recur 
          (vec (remove-trailing (dec dim) (take length (lazy-concat-repeat current))))
          (conj result current))
    ))))

(defn standard-factorization
  "Standard factorization of a Lyndon word `word`,
  i.e. word = uv with v the smallest nontrivial right-factor."
  [word]
  ; XXX this is brute force
  (let [sorted (sort
                 (fn [x y] (compare (second x) (second y)))
                 (map
                   (fn [i] [i (apply str (drop i word))])
                   (range 1 (count word))))
        break (first (first sorted))]
    [ (take break word) (drop break word) ]))

(defn bracketing
  [word]
  (if (= 1 (count word))
    (first word)
    (let [sf (standard-factorization word)]
      [(bracketing (first sf)) (bracketing (second sf))])))


(defn lie-bracket [a b]
  (lc/lc-apply-bilinear-function
    (fn [x y] (lc/lc-add (product x y)
                         (lc/lc-multiply -1 (product y x))))
    a b))

(defn bracket->lie-bracket [br]
  (if (sequential? br)
    (let [left (bracket->lie-bracket (first br))
          right (bracket->lie-bracket (second br))]
      (lie-bracket left right))
    (lc/lc-lift (sa/->ConcatWord [br]))))
    ;(bracket->lie-bracket-helper br)))


(defn spy [x]
  (println "spy => " x)
  x)

(defn lyndon-basis
  "The Lyndon basis in the letters 0,..,`dim`-1 of length <= order,
  lexicographically sorted."
  [dim order]
  (->> (lyndon-words dim order)
      (sort)
      (map (comp bracket->lie-bracket bracketing))))

;(defn dual-lyndon-basis
;   TODO
;  )

(defn lie-element
  " `m` a map, where (k,v) means the k'th Lyndon basis element is to be summed up with coefficient v."
  [m dim order]
  (let [basis (vec (lyndon-basis dim order))]
    (reduce
      lc/lc-add
      (map
        (fn [ [k v] ] (lc/lc-multiply v (basis k)))
        m))))

