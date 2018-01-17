(ns signature-invariants.invariants-test
  (:require 
    [clojure.test :refer :all]
    [signature-invariants.permutations :as perm]
    [signature-invariants.signature-test :as st]
    [hopf-algebra.shuffle-algebra :as sa]
    [hopf-algebra.linear-combination :as lc]
    [clatrix.core :as clatrix]
    [signature-invariants.signature :as signature]
    [hopf-algebra.hopf-algebra :refer [product to-str HopfAlgebra]]; coproduct antipode to-str to-latex]]
    [signature-invariants.invariants :refer :all]))

(defn- absolute-difference [x y]
  (if (or (clatrix/matrix? x)
          (clatrix/vec? x))
    (clatrix/norm (clatrix/- x y))
    (Math/abs (double (- x y)))))

(defn- close? [tolerance x y]
  (< (absolute-difference x y) tolerance))

(deftest matrix-stuff-test []
  (is (=
       2
       (clatrix/rank (clatrix/matrix [[1, 0, -1], [-2, 1, 4]]))))
  (is (=
       1
       (clatrix/rank (clatrix/matrix [[1, 0, -1], [2, 0, -2]]))))
  (is (=
       (clatrix/det (clatrix/matrix [ [1 2] [3 4] ]))
       (double (perm/det [ [1 2] [3 4] ]))))
  (is (close?
        0.00001
        (clatrix/det (clatrix/matrix [ [1 2 1] [4 5 6] [7 8 9] ] ))
        (double (perm/det [ [1 2 1] [4 5 6] [7 8 9] ])))))

(deftest misc
  )

(defn- sw [& args]
  (sa/->ShuffleWord args))
(deftest gl-invariants-test []
  (is
    (=
     [{(sw 0 1) 1, (sw 1 0) -1}]
     (gl-invariants 2 1)))

  (is
    (=
     [{(sw 0 1 0 1) 1,
       (sw 0 1 1 0) -1,
       (sw 1 0 0 1) -1
       (sw 1 0 1 0) 1}
      {(sw 0 0 1 1) 1,
       (sw 0 1 1 0) -1,
       (sw 1 0 0 1) -1
       (sw 1 1 0 0) 1}]
     (gl-invariants 2 2)))

  (testing "check dimensions"
    (doseq [weight (range 1 6)]
      (let [inv (gl-invariants 2 weight)]
        (is (=
             (/ (signature/n-choose-k (* 2 weight) weight) (inc weight)) ;- the Catalan number for weight
             (count inv)
             (count-independents inv 2 (* 2 weight))))))) ; the GL invariants are already linearly independent
  )

(deftest so-invariants-test []
  (testing "helper functions"
    (is
      (=
       ["c" "a" "b" "d"]
       (#'signature-invariants.invariants/multiply-specific-positions 4
                                    [ [[1 2] ["a" "b"]] [[0 3] ["c" "d"]] ])))

    (is
      (=
       { ["c" "a" "b" "d"] 10}
       (#'signature-invariants.invariants/lc-multiply-specific-positions 4
                                       [[1 2] [0 3]]
                                       [ { ["a" "b"] 2}
                                        { ["c" "d"] 5} ]))))

  (testing "dim=2"
    (doseq [k (range 1 4)]
      (is
        (=
         (* 2 (signature/n-choose-k (- (* 2 k) 1) (- k 1)))
         (count (so2-invariants (* 2 k)))
         (count-independents (so-invariants 2 (* 2 k)) 2 (* 2 k)))))))

(deftest permutation-invariants-test
  (is (=
       [{(sw 2) 1
         (sw 1) 1
         (sw 0) 1}]
       (permutation-invariants 3 1)))

  (is (=
       [{(sw 2 2) 1
         (sw 1 1) 1
         (sw 0 0) 1}
        {(sw 2 1) 1
         (sw 2 0) 1
         (sw 1 2) 1
         (sw 1 0) 1
         (sw 0 2) 1
         (sw 0 1) 1}]
       (permutation-invariants 3 2))))

;(deftest so2-invariants-jeremy-test
;
;  (is (=
;        '((0 0 1 1) (0 1 0 1) (0 1 1 0))
;        (possible-halves 4)))
;
;  )
