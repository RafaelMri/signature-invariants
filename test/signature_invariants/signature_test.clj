(ns signature-invariants.signature-test
  (:require [clojure.test :refer :all]
            [hopf-algebra.shuffle-algebra :as sa]
            [hopf-algebra.linear-combination :as lc]
            [signature-invariants.lyndon-words :as lw]
            [signature-invariants.signature :refer :all]))

(deftest dual-exponential-test []
  (testing "c"
    (is (= -1 (exp -1 1)))
    (is (= 1 (exp -1 2)))

    (is (= 2 (#'signature-invariants.signature/c 2 1)))
    (is (= -1/2 (#'signature-invariants.signature/c 2 2)))

    (is (= 3 (#'signature-invariants.signature/c 3 1)))
    (is (= -3/2 (#'signature-invariants.signature/c 3 2)))
    (is (= 1/3 (#'signature-invariants.signature/c 3 3))))

  (testing "de-shuffle"
    (is
      (=
       {(sa/->ShuffleWord [1 2]) 3,
        (sa/->ShuffleWord [2 1]) 1}
       (#'signature-invariants.signature/de-shuffle 2 (sa/->ShuffleWord [1 2])))))

  (is (=
       (dual-exponential (sa/->ShuffleWord [1 2]) 2)
       (dual-exponential (sa/->ShuffleWord [1 2]) 3)
       (dual-exponential (sa/->ShuffleWord [1 2]) 4)))
  )

(defn- cw [& w]
  (if (nil? w)
    (sa/->ConcatWord [])
    (sa/->ConcatWord w)))

(deftest exponential-test []
  (let [lie (lw/lie-element {0 1, 1 77, 2 555} 2 2)]
    (is
      (=
       { (cw) 1
         (cw 0) 1
         (cw 1) 77
         (cw 0 0) 1/2
         (cw 1 1) 5929/2
         (cw 0 1) 1187/2
         (cw 1 0) -1033/2}
       (exponential lie 2)))

    (is (= 555 (lc/lc-inner-product {(cw 0 1) 1/2, (cw 1 0) -1/2}
                                    (exponential lie 2))))
    ))

;(defn- inner-product-sw-cw [lc-1 lc-2]
  ;(lc/lc-inner-product (lw/sw->cw lc-1) lc-2))

;(deftest dual-basis-test []
;  (let [primal-basis (lw/primal-lyndon-basis 2 5)
;        dual-basis (lw/dual-lyndon-basis 2 5)
;        ; each Lie basis has its position as coefficient:
;        lie (lw/lie-element-cw (into {} (map (fn [x] [x (inc x)]) (range (count primal-basis)))))
;        group (exponential lie 5)
;        ]
;    (doseq [i (range 14)]
;      (is (= (inc i)
;             (inner-product-sw-cw
;               (dual-basis i)
;               lie)))
;      (is (= (inc i)
;             (inner-product-sw-cw
;               (lc/lc-apply-linear-function dual-exponential (dual-basis i))
;               group)))
;      )
;    )
;  )
