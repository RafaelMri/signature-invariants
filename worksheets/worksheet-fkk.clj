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
(ns worksheet-fkk
  (:require [gorilla-plot.core :as plot]))
(require '[signature-invariants.invariants :as inv] :reload)
(require '[hopf-algebra.linear-combination :as lc])
(require '[hopf-algebra.shuffle-algebra :as sa])
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(defn shuffle-concat [w-1 w-2]
  { (sa/->ShuffleWord (concat (:content w-1) (:content w-2))) 1})


(defn lc-shuffle-concat [lc-1 lc-2]
  (lc/lc-apply-bilinear-function
    shuffle-concat
    lc-1 lc-2))
  
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheet-fkk/lc-shuffle-concat</span>","value":"#'worksheet-fkk/lc-shuffle-concat"}
;; <=

;; @@
(def x { (sa/->ShuffleWord [0]) 1} )
(def y { (sa/->ShuffleWord [1]) 1} )

(def c lc-shuffle-concat)
(def s lc/lc-multiply)
(def *r lc/lc-multiply)



(def I1 (lc/lc-add
  (c x y)
  (*r -1/2 (s x y))))


(def I2 (lc/lc-add
  (s (c (s x y) y) x)
  (*r -1/2 (s (c (s x x) y) y))
  (*r -1/6 (s (s x x) (s y y)))))

(def I3 (lc/lc-add
	(s (c (s x y y) y) (s x x))
    (*r -1  (s (c (s x x y) y) (s x y)))
    (*r 1/3 (s (c (s x x x) y) (s y y)))
    (*r -1/12 (s x x x y y y))))


(lc/lc-to-str I1)
(lc/lc-to-str I2)
(lc/lc-to-str I3)


(def d2w2 (inv/gl-invariants 2 2))
(def d2w3 (inv/gl-invariants 2 3))

;(type d2w2)

(assert (=
(inv/count-independents d2w2 2 4)
(inv/count-independents (conj d2w2 I2) 2 4)))

(assert (=
(inv/count-independents d2w3 2 6)
(inv/count-independents (conj d2w3 I3) 2 6)))

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@

;; @@

;; @@

;; @@
