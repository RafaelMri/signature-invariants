(defproject signature-invariants "0.1.0-SNAPSHOT"
  :description "A library to calculate the invariants of `Diehl/Reizenstein - Invariants of multidimensional signals based on their signature, 2018.`"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [hopf-algebra "0.1.0-SNAPSHOT"]
                 [org.clojure/math.combinatorics "0.1.4"]
                 [clatrix "0.5.0"]
                 [spyscope "0.1.5"]
                 ]
  :target-path "target/%s"
  :plugins
    [[lein-gorilla "0.4.0"]]
  :profiles {:uberjar {:aot :all}
             :repl {:dependencies [[gorilla-renderable "1.0.0"]
                                   [org.clojure/tools.namespace "0.2.11"]]}
             :test {:dependencies [[gorilla-renderable "1.0.0"]]}})
