(ns signature-invariants.lyndon-words-test
  (:require 
    [clojure.test :refer :all]
    [hopf-algebra.shuffle-algebra :as sa]
    [hopf-algebra.linear-combination :as lc]
    [hopf-algebra.linear-combination-plus :as lcp]
    [signature-invariants.lyndon-words :refer :all]))


; Lyndon basis
(def x1 {(sa/->ConcatWord [0]) 1})
(def x2 {(sa/->ConcatWord [1]) 1})
(def b lie-bracket)

(def first-5-levels
  [
  ; 0 -> 1
   x1
  ; 1 -> 2
   x2
  ; 2 -> [1,2]
   (b x1 x2)
  ; 3 -> [[1,2],2]
   (b (b x1 x2) x2)
  ; 4 -> [1,[1,2]
   (b x1 (b x1 x2))
  ; 5 -> [[[1,2],2],2]
   (b (b (b x1 x2) x2) x2)
  ; 6 -> [1,[[1,2],2]]
   (b x1 (b (b x1 x2) x2))
  ; 7 -> [1,[1,[1,2]]]
   (b x1 (b x1 (b x1 x2)))
  ; 8 -> [ [[[1,2],2],2], 2]
   (b (b (b (b x1 x2) x2) x2) x2)
  ; 9 -> [[1,[1,2], [1,2]
   (b (b x1 (b x1 x2)) (b x1 x2))
  ; 10 -> [[1,2], [[1,2],2]]
   (b (b x1 x2) (b (b x1 x2) x2))
  ; 11 -> [1, [[[1,2],2],2]]
   (b x1 (b (b (b x1 x2) x2) x2))
  ; 12 -> [1, [1,[[1,2],2]]]
   (b x1 (b x1 (b (b x1 x2) x2)))
  ; 13 -> [1, [1,[1,[1,2]]]]
   (b x1 (b x1 (b x1 (b x1 x2))))])



(deftest misc-test []
  (is (=
       '(1 2 3 1 2 3 1 2 3 1)
       (take 10 (lazy-concat-repeat [1 2 3]))))

  (is (=
       '( 1 4 2 3 )
       (remove-trailing 4 '( 1 4 2 3 4 4 ))))
  (is (=
       '( )
       (remove-trailing 4 '( 4 4 ))))

  (is (=
       #{ [0] [0, 0, 1] [0, 1] [0, 1, 1] [1] }
       (into #{} (lyndon-words 2 3))))

  (is (=
       #{ [0] [0, 0, 0, 1] [0, 0, 1] [0, 0, 1, 1] [0, 1] [0, 1, 1] [0, 1, 1, 1] [1] }
       (into #{} (lyndon-words 2 4))))

  (is (=
      [ [1] [1 2] ]
      (standard-factorization [1 1 2])))

  (is (=
      [ [1] [ 2] ]
      (standard-factorization [1 2])))

  (is (=
      [1 [1 2]]
      (bracketing [1 1 2])))

  (doseq [ [word expected-bracket]
          {
          [1 2] [1, 2]
          [1 1 2] [1, [1, 2]]
          [1 2 2] [[1, 2], 2]
          [1 1 1 2] [1, [1, [1, 2]]]
          [1 1 2 2] [1, [[1, 2], 2]]
          [1 2 2 2] [[[1, 2], 2], 2]
          [1 1 1 1 2] [1, [1, [1, [1, 2]]]]
          [1 1 1 2 2] [1, [1, [[1, 2], 2]]]
          [1 1 2 1 2] [[1, [1, 2]], [1, 2]]
          [1 1 2 2 2] [1, [[[1, 2], 2], 2]]
          [1 2 1 2 2] [[1, 2], [[1, 2], 2]]
          [1 2 2 2 2] [[[[1, 2], 2], 2], 2] } ]
    (testing (str "word " word)
      (is (=
           expected-bracket
           (bracketing word)))
      )
    )

  (is (= {(sa/->ConcatWord [0 1]) 1, (sa/->ConcatWord [1 0]) -1}
         (bracket->lie-bracket [0 1])))

  (is (= (into #{} first-5-levels)
         (into #{} (lyndon-basis 2 5))))
  )

(deftest chen-fox-lyndon-breakpoints-test
  (is (=
       [ [1, 3, 3, 2] [1, 1, 3] [1] ] 
       (chen-fox-lyndon-factorization [1,3,3,2,1,1,3,1])))

  (let [words (sort (lyndon-words 2 7)),
        lie
          (reduce
            lcp/+
            (map-indexed
              (fn [i w]
                (-> w
                    (bracketing)
                    (bracket->lie-bracket)
                    (lcp/lc->lcp (lcp/->FreeCommutative { (str "a" i) 1}))))
                words))]
    (doseq [i (range (count words))]
      (is (=
           { (lcp/->FreeCommutative { (str "a" i) 1}) 1 }
           (lcp/lc-lc-plus-inner-product (sa/sw->cw (word->dual-basis (nth words i))) lie)))
      )))
