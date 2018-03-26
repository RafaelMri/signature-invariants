(ns signature-invariants.signature-test
  (:require [clojure.test :refer :all]
            [hopf-algebra.shuffle-algebra :as sa]
            [hopf-algebra.linear-combination :as lc]
            [hopf-algebra.linear-combination-plus :as lcp]
            [signature-invariants.lyndon-words :as lw]
            [signature-invariants.signature :refer :all]))


(deftest dual-logarithm-test
  (is (=
    [ [[:a :b :c :d]] [[:a] [:b :c :d]] [[:a :b] [:c :d]] [[:a :b :c] [:d]]
      [[:a] [:b] [:c :d]] [[:a] [:b :c] [:d]] [[:a :b] [:c] [:d]] [[:a] [:b] [:c] [:d]] ]
    (#'signature-invariants.signature/all-splits [:a :b :c :d])))

  (let [words (sort (lw/lyndon-words 2 4)),
        lie
          (reduce
            lcp/+
            (map-indexed
              (fn [i w]
                (-> w
                    (lw/word->basis)
                    (lcp/lc->lcp (lcp/->FreeCommutative { (str "a" i) 1}))))
                words))
        exp (time (exponential-symbolic lie {(lcp/->FreeCommutative {}) 1}))] ; XXX very slow
    (doseq [i (range (count words))]
      (is (=
           { (lcp/->FreeCommutative { (str "a" i) 1}) 1 }
           (lcp/lc-lc-plus-inner-product (sa/sw->cw (lw/word->dual-basis (nth words i))) lie)
           (lcp/lc-lc-plus-inner-product
              (sa/sw->cw (lc-dual-logarithm (lw/word->dual-basis (nth words i))))
              exp))))))

;(defn- cw [& w]
;  (if (nil? w)
;    (sa/->ConcatWord [])
;    (sa/->ConcatWord w)))

;(deftest exponential-test []
;  (let [lie (lw/lie-element {0 1, 1 77, 2 555} 2 2)]
;    (is
;      (=
;       { (cw) 1
;         (cw 0) 1
;         (cw 1) 77
;         (cw 0 0) 1/2
;         (cw 1 1) 5929/2
;         (cw 0 1) 1187/2
;         (cw 1 0) -1033/2}
;       (exponential lie 2)))
;
;    (is (= 555 (lc/inner-product {(cw 0 1) 1/2, (cw 1 0) -1/2}
;                                    (exponential lie 2)))))
;
;
;  (let [fa (fn [i] {(lcp/->FreeCommutative {(str "a" i) 1}) 1})
;        lie (lw/lie-element-symbolic {0 (fa 0), 1 (fa 1), 2 (fa 2)} 2 2)
;        unit {(lcp/->FreeCommutative {}) 1}]
;    (is (=
;         (fa 2)
;         (lcp/lc-lc-plus-inner-product 
;           {(cw 0 1) 1/2, (cw 1 0) -1/2}
;           (exponential-symbolic lie unit))))
;      )
;
;  (let [lie {(sa/->ConcatWord [1]) {(lcp/->FreeCommutative {"a" 1}) 1}}
;        unit {(lcp/->FreeCommutative {}) 1}]
;    (is (=
;         {(sa/->ConcatWord []) unit,
;          (sa/->ConcatWord [1]) {(lcp/->FreeCommutative {"a" 1}) 1}
;          (sa/->ConcatWord [1 1]) {(lcp/->FreeCommutative {"a" 2}) 1/2}}
;        (exponential-symbolic lie unit 2))))
;  )
