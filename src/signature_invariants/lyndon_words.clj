(ns
  signature-invariants.lyndon-words
  "Methods to calculate Lyndon words and the Lyndon basis of the free Lie algebra."
  (:require [signature-invariants.young-tableau :as yt]
            [hopf-algebra.linear-combination :as lc]
            [hopf-algebra.linear-combination-plus :as lcp]
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
          (conj result current))))))

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
  (lc/apply-bilinear-function
    (fn [x y] (lc/+ (product x y)
                         (lc/* -1 (product y x))))
    a b))

(defn bracket->lie-bracket
  "[1 2] => 12 - 21"
  [br]
  (if (sequential? br)
    (let [left (bracket->lie-bracket (first br))
          right (bracket->lie-bracket (second br))]
      (lie-bracket left right))
    (lc/lift (sa/->ConcatWord [br]))))
    ;(bracket->lie-bracket-helper br)))

(defn word->basis
  "Lyndon word to corresponding Lie basis element."
  [word]
  (-> word
      (bracketing)
      (bracket->lie-bracket)))

(defn lyndon-basis
  "The Lyndon basis in the letters 0,..,`dim`-1 of length <= order,
  lexicographically sorted."
  [dim order]
  (->> (lyndon-words dim order)
      (sort)
      (map word->basis)))

(defn- chen-fox-lyndon-breakpoints-next-ij [k word]
  (loop [i k,
         j (inc k)]
    (if (or
          (>= j (count word))
          (> (nth word i) (nth word j)))
      [i j]
      (recur
        (if (= (nth word i) (nth word j)) (inc i) k)
        (inc j)))))

(defn- chen-fox-lyndon-breakpoints-next-ks [k word]
  (let [ [i j] (chen-fox-lyndon-breakpoints-next-ij k word)]
    (loop [k' (+ k (- j i)), 
           result [k']]
      (if (>= k' (inc i))
        result
        (recur
          (+ k' (- j i))
          (conj result k'))))))

(defn chen-fox-lyndon-breakpoints
  "The algorithm from Duval - 1983 - Factorizing words over and ordered alphabet."
  [word]
  (loop [k 0,
         result [] ]
    (if (>= k (count word))
      result
      (let [ks (chen-fox-lyndon-breakpoints-next-ks k word)]
        (recur
          (last ks)
          (concat result ks))))))

(defn- subword [coll n m]
  (take
    (- m n)
    (drop n coll)))

(defn chen-fox-lyndon-factorization [word]
  (loop [i 0,
         breakpoints (chen-fox-lyndon-breakpoints word)
         result []
         ]
    (if (empty? breakpoints)
      result
      (recur
        (first breakpoints)
        (rest breakpoints)
        (conj result
              (subword word i (first breakpoints)))))))

(defn- fact
  "n !"
  [n]
  (apply * (range 1 (inc n))))

(defn word->dual-basis
  "Word to corresponding dual PWB element; Reutenauer - Thm 5.3.
   In particular for a Lyndon word, gives the dual to the corresponding Lie element."
  [word]
  (if (empty? word)
    {(sa/->ShuffleWord []) 1}
    (let [cfl (chen-fox-lyndon-factorization word)]
    (if (= 1 (count cfl)) ; Lyndon
      (lc/apply-linear-function
        (fn [sw] {(sa/->ShuffleWord (concat [(first word)] (:content sw))) 1})
        (word->dual-basis (rest word)))
      (lc/*
        (/ 1 (apply * (map (fn [[_ f]] (fact f)) (frequencies cfl))))
        (apply lc/* (map word->dual-basis cfl)))))))

(defn lie-element
  " `m` a map, where (k,v) means the k'th Lyndon basis element is to be summed up with coefficient v."
  [m dim order]
  (let [basis (vec (lyndon-basis dim order))]
    (reduce
      lc/+
      (map
        (fn [ [k v] ] (lc/* v (basis k)))
        m))))

(defn lie-element-symbolic
  " `m` a map, where (k,v) means the k'th Lyndon basis element is to be summed up with coefficient v."
  [m dim order]
  (let [basis (vec (lyndon-basis dim order))]
    (reduce
      lcp/+
      (map
        (fn [ [k v] ]
          (lcp/apply-linear-function'
            (fn [x] {x v})
            (basis k)))
        m))))

