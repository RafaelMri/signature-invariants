;; gorilla-repl.fileformat = 1

;; **
;;; # Gorilla REPL
;;; 
;;; Welcome to gorilla :-)
;;; 
;;; Shift + enter evaluates code. Hit ctrl+g twice in quick succession or click the menu icon (upper-right corner) for more commands ...
;;; 
;;; It's a good habit to run each worksheet in its own namespace: feel free to use the declaration we've provided below if you'd like.
;; **

;; @@
(ns talented-volcano
  (:require [gorilla-plot.core :as plot]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(require '[signature-invariants.invariants :as inv] :reload)
(require '[hopf-algebra.shuffle-algebra :as sa])

(require '[hopf-algebra.linear-combination :as lc] :reload)

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(lc/lc-to-str (first (inv/gl-invariants 4 1)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-string'>&quot;- 3012 - 1023 + 3210 + 2301 - 2310 + 0312 - 3201 - 0132 - 1302 + 1032 - 2103 + 0123 + 2013 - 0321 + 0231 + 1320 - 2031 + 1203 + 3102 - 0213 + 2130 - 3120 - 1230 + 3021&quot;</span>","value":"\"- 3012 - 1023 + 3210 + 2301 - 2310 + 0312 - 3201 - 0132 - 1302 + 1032 - 2103 + 0123 + 2013 - 0321 + 0231 + 1320 - 2031 + 1203 + 3102 - 0213 + 2130 - 3120 - 1230 + 3021\""}
;; <=

;; @@
(def P (for [i (range 4)]
  (for [j (range 4)]
    (if (= i j)
      {}
       { (sa/->ShuffleWord [i j]) 1, (sa/->ShuffleWord [j i]) -1})
    )))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;talented-volcano/P</span>","value":"#'talented-volcano/P"}
;; <=

;; @@
(defn p [i j]
  (nth (nth P i) j))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;talented-volcano/p</span>","value":"#'talented-volcano/p"}
;; <=

;; @@
(p 0 2)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"list-like","open":"<span class='clj-record'>#hopf_algebra.shuffle_algebra.ShuffleWord{</span>","close":"<span class='clj-record'>}</span>","separator":" ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:content</span>","value":":content"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[0 2]"}],"value":"[:content [0 2]]"}],"value":"#hopf_algebra.shuffle_algebra.ShuffleWord{:content [0 2]}"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[#hopf_algebra.shuffle_algebra.ShuffleWord{:content [0 2]} 1]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"list-like","open":"<span class='clj-record'>#hopf_algebra.shuffle_algebra.ShuffleWord{</span>","close":"<span class='clj-record'>}</span>","separator":" ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:content</span>","value":":content"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[2 0]"}],"value":"[:content [2 0]]"}],"value":"#hopf_algebra.shuffle_algebra.ShuffleWord{:content [2 0]}"},{"type":"html","content":"<span class='clj-long'>-1</span>","value":"-1"}],"value":"[#hopf_algebra.shuffle_algebra.ShuffleWord{:content [2 0]} -1]"}],"value":"{#hopf_algebra.shuffle_algebra.ShuffleWord{:content [0 2]} 1, #hopf_algebra.shuffle_algebra.ShuffleWord{:content [2 0]} -1}"}
;; <=

;; @@
; (4) in Luque, Thibon - 2002 - Pfaffian and Hafnian identities in shuffle algebras
(lc/lc-add
  (lc/lc-multiply (p 0 1) (p 2 3))
  (lc/lc-multiply -1 (p 0 2) (p 1 3))
  (lc/lc-multiply (p 1 2) (p 0 3))
  (lc/lc-multiply -1 (first (inv/gl-invariants 4 1)))
  )
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[],"value":"{}"}
;; <=

;; @@

;; @@
