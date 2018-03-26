(ns
  signature-invariants.worksheets.worksheet-invariants
  (:require [hopf-algebra.shuffle-algebra :as sa]
            [hopf-algebra.linear-combination :as lc]
            [hopf-algebra.linear-combination-plus :as lcp]
            [hopf-algebra.hopf-algebra :refer [product]]; coproduct antipode to-str to-latex]]
            [clojure.math.numeric-tower :refer [expt]]
            [signature-invariants.lyndon-words :as lw]
            [spyscope.core]
            ))

(defn print-lyndon-basis [dim order]

    (doseq [word (sort (lw/lyndon-words dim order))]
      (println word)
      (println "bracketing=" (lw/bracketing word))
      (println "lie element=" (lc/to-str (lw/word->basis word)))
      (println "dual lie element=" (lc/to-str (lw/word->dual-basis word))))

  )




(defn run []
  ;(print-lyndon-basis 2 5)
  ;..

  )
