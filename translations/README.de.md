# mrwin: Kausale Win-Statistiken für hierarchische zusammengesetzte Endpunkte

[English](../README.md) | [Tiếng Việt](README.vi.md) | [中文](README.zh.md) | [日本語](README.ja.md) | [한국어](README.ko.md) | [Français](README.fr.md) | [Deutsch](README.de.md)

## Was dieses Paket tut

`mrwin` ist ein R-Paket, das eine spezifische Frage in der Epidemiologie beantwortet:

**Führt eine genetisch vorhergesagte Exposition dazu, dass Menschen bessere oder
schlechtere Gesamtergebnisprofile aufweisen, wenn die Ergebnisse klinisch nach
Schweregrad geordnet sind?**

Bei vielen Erkrankungen erleben Patienten nicht nur ein einzelnes Ereignis. Eine
Person mit kardiorenaler Erkrankung kann sterben, wegen Herzinsuffizienz
hospitalisiert werden oder eine Nierenfunktionsverschlechterung entwickeln.
Diese Ereignisse sind nicht austauschbar: Tod ist schlimmer als
Krankenhauseinweisung, was schlimmer ist als eine Laboranomalie.
Standardanalysen behandeln alle Ereignisse gleich oder analysieren sie einzeln,
wobei die klinische Rangfolge verloren geht, die für Patienten und Kliniker am
wichtigsten ist.

`mrwin` kombiniert zwei statistische Rahmenwerke, um dies zu adressieren:

1. **Win-Statistiken** -- paarweise Vergleiche, die klinische
   Schweregradhierarchien respektieren (Tod schlägt Krankenhauseinweisung
   schlägt Biomarker-Verschlechterung).
2. **Mendelsche Randomisierung** -- genetische instrumentelle Variablen, die
   kausale Aussagen unterstützen, selbst wenn Confounder nicht gemessen werden.

Das Ergebnis ist eine **dosisstandardisierte kausale Win-Ratio (DS-CWR)**: eine
einzige Zahl, die zusammenfasst, ob eine genetisch vorhergesagte Exposition das
priorisierte Ergebnisprofil besser oder schlechter macht.

---

## Mendelsche Randomisierung: Eine kurze Einführung

### Das Problem: Nicht gemessenes Confounding

In der beobachtenden Epidemiologie möchten wir oft wissen, ob eine Exposition
(z. B. LDL-Cholesterin) ein Ergebnis (z. B. Herzinsuffizienz) verursacht. Die
Schwierigkeit besteht darin, dass Menschen mit hohem LDL sich auch in vielerlei
anderer Hinsicht von Menschen mit niedrigem LDL unterscheiden -- Ernährung,
Bewegung, sozioökonomischer Status, andere Medikamente. Dies sind Confounder.
Selbst nach Adjustierung für gemessene Confounder können nicht gemessene
Confounder verbleiben.

### Die Lösung: Genetik als natürliches Experiment

Bei der Empfängnis erbt jede Person zufällig ein Allel von jedem Elternteil an
jedem genetischen Locus. Dies wird **Mendelsche Segregation** genannt. Wenn eine
genetische Variante eine Exposition beeinflusst (z. B. eine Variante im
LDLR-Gen erhöht das LDL-Cholesterin), dann sind Personen, die diese Variante
tragen, im Durchschnitt lebenslang höherem LDL ausgesetzt -- nicht wegen ihrer
Ernährung oder ihres Lebensstils, sondern wegen ihres Genotyps.

Diese zufällige Zuteilung ist die Grundlage der **Mendelschen Randomisierung
(MR)**. Durch die Verwendung genetischer Varianten als **instrumentelle
Variablen** schätzt MR den kausalen Effekt der Exposition auf das Ergebnis und
umgeht dabei nicht gemessenes Confounding.

### Kernannahmen der MR

Damit das genetische Instrument gültig ist, müssen drei Bedingungen erfüllt sein:

1. **Relevanz**: Die genetische Variante muss mit der Exposition assoziiert sein.
2. **Unabhängigkeit**: Die genetische Variante darf nicht mit Confoundern
   assoziiert sein (dies folgt aus der Mendelschen Segregation in gut
   konzipierten Studien).
3. **Exklusionsrestriktion**: Die genetische Variante darf das Ergebnis nur
   über die Exposition beeinflussen, nicht über andere Pfade.

### Datentypen in der MR

| Datentyp | Was es ist | Beispiel |
|---|---|---|
| Individualebene | Eine Zeile pro Person mit Genotyp, Exposition, Ergebnis | Eine Kohorte von 10.000 Personen mit SNP-Daten, LDL-Messungen und Krankenhausakten |
| Summary-Level | Pro-SNP-Effektschätzungen aus GWAS | „SNP rs12345 hat Effekt 0,05 (SE 0,01) auf LDL" |
| Polygenetischer Risikoscore (PRS) | Gewichtete Summe vieler SNPs | PRS für LDL-Cholesterin, berechnet aus 100 SNPs |

`mrwin` verwendet **Daten auf Individualebene** mit einem **PRS-Instrument**.
Der PRS dient als kontinuierliches Instrument, das die kumulative genetische
Prädisposition für die Exposition erfasst.

---

## Win-Statistiken: Der Ziel-Estimand

### Was ist eine Win-Statistik?

Eine Win-Statistik vergleicht jedes Paar von Personen in einer Studie. Für jedes
Paar (A, B) folgt der Vergleich einer klinischen Hierarchie:

1. Zuerst Vergleich des wichtigsten Ergebnisses (z. B. Tod). Wenn A gestorben
   ist und B nicht, gewinnt B.
2. Bei Gleichstand beim ersten Ergebnis, gehe zum nächsten über (z. B.
   Krankenhauseinweisung).
3. Fahre in der Hierarchie fort, bis eine Person gewinnt oder alle Ergebnisse
   gleich sind.

Die **Win-Ratio** ist die Anzahl der Gewinne für die höher exponierte Gruppe
geteilt durch die Anzahl der Gewinne für die niedriger exponierte Gruppe.

### Warum Win-Statistiken wichtig sind

| Methode | Behandelt Ereignisse gleich? | Respektiert klinische Rangfolge? | Berücksichtigt konkurrierende Risiken? |
|---|---|---|---|
| Time-to-first-event | Ja | Nein | Nein |
| Analyse zusammengesetzter Endpunkte | Ja | Nein | Teilweise |
| Win-Ratio / Win-Statistik | Nein | Ja | Ja |

Win-Statistiken wurden von Pocock et al. (2012) eingeführt und wurden in über 36
randomisierten klinischen Studien zwischen 2022 und 2024 eingesetzt.

### Die kausale Win-Ratio (cCWR)

Standard-Win-Statistiken erfordern Randomisierung oder starke Ignorierbarkeit.
In der beobachtenden Epidemiologie ist diese Annahme selten vertretbar.

`mrwin` definiert die **kontinuierliche kausale Win-Ratio (cCWR)** als
Gradienten-Estimand über Quantile des genetischen Instruments:

- Teile die Population in geordnete Strata nach PRS (z. B. Dezile).
- Innerhalb jedes benachbarten Stratumpaares berechne die Win-Ratio.
- Poole über Strata mittels instrumentenstandardisierter Gradienten und
  GLS-Metaanalyse.

Das Ergebnis ist ein **marginaler, populationsweiter kausaler Effekt**, der
keine starke Ignorierbarkeit erfordert.

---

## Kernformeln: Erklärung in einfacher Sprache

### 1. Der hierarchische Vergleichskernel

Für zwei Personen *i* und *j*, vergleiche Ergebnisse von höchster zu niedrigster
Priorität:

- Wenn *j* das Ereignis bei Priorität *k* hatte, aber *i* nicht (oder *i*
  länger überlebt hat): *i* gewinnt (+1).
- Wenn *i* das Ereignis bei Priorität *k* hatte, aber *j* nicht (oder *j*
  länger überlebt hat): *i* verliert (-1).
- Bei Gleichstand: gehe zur nächsten Priorität.

**In Worten**: „Hatte diese Person ein besseres Ergebnisprofil als die andere,
unter Berücksichtigung der klinischen Schweregrad-Reihenfolge?"

### 2. Die stratumspezifische Win-Ratio

Innerhalb eines Stratums (z. B. das 6. PRS-Dezil vs. das 5.):

```
theta_d = (Anzahl der Gewinne für das höhere Stratum) / (Anzahl der Verluste)
log_theta_d = log(theta_d)
```

**In Worten**: „Haben Personen mit etwas höherer genetischer Prädisposition für
die Exposition tendenziell häufiger gewonnen als verloren?"

### 3. Der instrumentenstandardisierte Gradient (ISG)

```
Delta_X_d = Mittelwert(Exposition in Stratum d) - Mittelwert(Exposition in Stratum d-1)
ISG_d = log_theta_d / Delta_X_d
```

**In Worten**: „Wie stark ändert sich die Log-Win-Ratio pro Einheitszunahme der
Exposition, wie durch Genetik vorhergesagt?"

### 4. GLS-Pooling

Die ISG-Werte über alle Stratumpaare werden mittels Generalisierter
Kleinste-Quadrate-Schätzung mit einem geschrumpften Kovarianzschätzer gepoolt:

```
delta_GLS = gepoolter ISG über alle Strata
DS-CWR = exp(delta_GLS)
```

**In Worten**: „Was ist der gesamte kausale Effekt der Exposition auf das
priorisierte Ergebnisprofil?"

- **DS-CWR > 1**: Die Exposition ist mit einem besseren Ergebnisprofil
  assoziiert.
- **DS-CWR < 1**: Die Exposition ist mit einem schlechteren Ergebnisprofil
  assoziiert.
- **DS-CWR = 1**: Kein Hinweis auf einen kausalen Effekt.

### 5. Der Multiplikator-Bootstrap

Unsicherheit wird geschätzt durch:
1. Perturbation der GWAS-Gewichte (externe Unsicherheit).
2. Ziehen zufälliger Multiplikator-Gewichte (interne Unsicherheit).
3. Neuberechnung der gesamten Pipeline für jede Bootstrap-Iteration.

Dies ergibt ein Konfidenzintervall, das beide Unsicherheitsquellen
berücksichtigt.

### 6. Die Step-Down-Pleiotropie-Diagnostik (SDPD)

Die SDPD testet, ob das genetische Instrument das Ergebnis mit der höchsten
Priorität über andere Pfade als die Exposition beeinflusst (direkte
Pleiotropie). Dies geschieht durch Anpassung einer MR-Egger-Regression auf die
per-SNP Expositions- und Ergebniseffekte.

**In Worten**: „Ist das genetische Instrument valide, oder beeinflusst es das
Ergebnis über andere Pfade?"

---

## Klinische Beispiele

### Beispiel 1: Kardiovaskuläre Erkrankung

**Forschungsfrage**: Verschlechtert genetisch vorhergesagtes LDL-Cholesterin den
gesamten kardiorenalen Verlauf?

**Endpunkthierarchie** (höchste zu niedrigste Priorität):
1. Tod
2. Krankenhauseinweisung wegen Herzinsuffizienz
3. Nierenfunktionsverschlechterung

**Interpretation**: Wenn DS-CWR = 0,82 (95%-KI: 0,70 bis 0,96), bedeutet dies,
dass pro Einheitszunahme des genetisch vorhergesagten LDL das priorisierte
kardiorenale Ergebnisprofil um 18 % schlechter wird. Das KI kreuzt nicht 1, was
auf einen statistisch signifikanten schädlichen Effekt hindeutet.

**Klinische Bedeutung**: Höheres LDL verursacht nicht nur Herzinfarkte -- es
verschiebt den gesamten Krankheitsverlauf in Richtung früheren Todes, mehr
Krankenhauseinweisungen und schnelleren Nierenfunktionsverlust, in dieser
klinischen Reihenfolge.

### Beispiel 2: Demenz

**Forschungsfrage**: Führt genetisch vorhergesagte Schlafstörung dazu, dass
Menschen früher schwere neurokognitive Folgen erleiden?

**Endpunkthierarchie**:
1. Tod
2. Demenzdiagnose
3. Pflegeheimeinweisung

**Interpretation**: Wenn DS-CWR = 0,91 (95%-KI: 0,78 bis 1,06), deutet die
Punktschätzung auf Schaden hin, aber das KI kreuzt 1. Es gibt keine starke
Evidenz für einen kausalen Effekt.

**Warum das wichtig ist**: In der Standard-MR von Schlaf auf Demenz werden
Personen, die sterben, bevor sie Demenz entwickeln, aus der Analyse
ausgeschlossen. Dies ist eine Form von Survivor Bias: Eine schädliche
Schlafvariante kann protektiv erscheinen, weil ihre Träger nie lange genug
leben, um diagnostiziert zu werden. `mrwin` vermeidet dies, indem Tod als das
Ereignis mit der höchsten Priorität behandelt wird, sodass differentielles
Überleben in den Estimanden aufgenommen und nicht wegkonditioniert wird.

### Beispiel 3: Krebsepidemiologie

**Forschungsfrage**: Verschlechtert genetisch vorhergesagter BMI die
priorisierten Krebsergebnisse?

**Endpunkthierarchie**:
1. Krebstod
2. Progression oder Metastasierung
3. Rezidiv

**Interpretation**: Wenn DS-CWR = 0,75 (95%-KI: 0,60 bis 0,94), bedeutet dies,
dass ein höherer genetisch vorhergesagter BMI den Krebsverlauf auf jeder Ebene
der Hierarchie in Richtung schlechterer Ergebnisse verschiebt.

**Klinische Bedeutung**: BMI erhöht nicht nur das Krebsrisiko -- er
verschlechtert den gesamten Krankheitsverlauf vom Rezidiv bis zum Tod.

---

## Installation

```r
# Von GitHub (Entwicklungsversion)
# install.packages("devtools")
devtools::install_github("nguyenminh2301/mrwin")
```

---

## Schnellstart

### Simulierte Daten (keine echte Kohorte erforderlich)

```r
library(mrwin)

# Simuliere einen kardiorenalen Datensatz
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)

# Passe das Modell an
fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)

# Ergebnisse anzeigen
print(fit)
summary(fit)
plot(fit)
```

### Die Ausgabe lesen

```
mrwin fit
  N: 500
  SNPs: 20
  Priorities: 3
  Strata: 5
  Bootstrap: 100 (100 valid)
  Backend: dense
  Adjustment: none
  delta_GLS: -0.2000
  DS-CWR: 0.8187
  95% CI: 0.7000 to 0.9600
  Q: 1.2345 ( df = 2 , p = 0.5394 )
```

| Feld | Was es bedeutet | Worauf zu achten ist |
|---|---|---|
| `delta_GLS` | Gepoolter Effekt auf der Log-Skala | Negativ = schädlich, positiv = vorteilhaft |
| `DS-CWR` | Exponentierter Effekt (Win-Ratio) | < 1 = schädlich, > 1 = vorteilhaft, 1 = kein Effekt |
| `95% CI` | Bootstrap-Konfidenzintervall | Kreuzt es 1? |
| `Q` | Heterogenitätsstatistik | Kleiner p-Wert = Effekt variiert über Strata |
| `Warnings` | Strukturierte Vorbehalte | Immer vor der Interpretation prüfen |

### Das Ergebnis interpretieren

Wenn DS-CWR = 0,82 (95%-KI: 0,70 bis 0,96):

- Die Exposition ist mit einem um 18 % schlechteren priorisierten
  Ergebnisprofil assoziiert.
- Das KI kreuzt nicht 1, der Effekt ist also statistisch signifikant.
- Prüfen Sie den Q-p-Wert: Wenn signifikant, kann der Effekt über den
  Expositionsbereich variieren.
- Prüfen Sie die Warnungen: Wenn `sdpd_rejected`, kann Pleiotropie die
  Schätzung invalidieren.

---

## Vorlage für reale Kohorten

```r
snp_cols <- paste0("rs", 1:20)

fit <- mrwin(
  data = dat,
  endpoint = mrwin_endpoint(
    time = c("t_death", "t_hf", "t_renal"),
    status = c("d_death", "d_hf", "d_renal"),
    priority = c("death", "heart_failure", "renal_decline")
  ),
  genotype = snp_cols,
  exposure = "ldl_cholesterol",
  gwas = mrwin_gwas(beta = beta_hat, se = se_beta, snp = snp_cols),
  covariates = c("age", "sex", "PC1", "PC2"),
  controls = mrwin_controls(
    n_strata = 10,
    bootstrap = 500,
    seed = 20260510,
    adjustment = "ordinal_iptw",
    sdpd_scale = "both"
  )
)

summary(fit)
plot(fit, type = "forest")
mrwin_report(fit, file = "analysis_report.md", format = "markdown")
tidy(fit)
```

---

## Ausgabeformate

### Print

```r
print(fit)
```

Zeigt das Hauptergebnis: Stichprobengröße, SNPs, Strata, Bootstrap-Validität,
DS-CWR, KI, Q-Statistik und Warnungen.

### Summary

```r
summary(fit)
```

Erstellt einen formatierten Bericht mit:
- DS-CWR-Hauptschätzung, primärem (Fieller-)KI, p-Wert und optionalem
  Pleiotropie-begrenztem KI
- Benachbarten Stratumgradienten (log-theta, Delta-X, ISG, CWR)
- Heterogenitäts-Q-Statistik
- Dem Delta-Methoden-Intervall als beschriftete Referenz (das Fieller-Intervall
  ist das primäre; siehe *Konfidenzintervalle* unten)
- SDPD-Diagnostik (wenn aktiviert)
- Strukturierte Warnungen

### Tidy (für Tabellen und weitere Analysen)

```r
tidy(fit)
#   term estimate    delta se_delta statistic p_value  ci_low ci_high ...
# DS-CWR    0.82 -0.2000   0.075    -2.667  0.0077  0.7000  0.9600 ...
```

### Report

```r
mrwin_report(fit, format = "markdown")
mrwin_report(fit, file = "report.md", format = "markdown")
```

### Plots

```r
plot(fit, type = "isg")       # Benachbarte ISG-Gradienten (Standard)
plot(fit, type = "forest")    # Forest-Plot der benachbarten CWR
plot(fit, type = "bootstrap") # Bootstrap-Verteilung von log-theta
```

---

## Kovariaten-Adjustierung

Wenn Confounder gemessen werden (Alter, Geschlecht, Hauptkomponenten),
verwenden Sie ordinale IPTW-Adjustierung:

```r
fit_adj <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  covariates = cbind(u_proxy = scale(dat$U)),
  controls = mrwin_controls(
    n_strata = 5,
    bootstrap = 100,
    seed = 2,
    adjustment = "ordinal_iptw",
    sdpd_scale = "both",
    pleiotropy_bias_radius = 0.05
  )
)
```

Die IPTW-Adjustierung verbessert die Präzision, indem Kovariaten über
PRS-Strata balanciert werden. Die SDPD testet auf direkte Pleiotropie.

---

## Stärken und Limitationen

### Stärken

- **Kausale Inferenz bei nicht gemessenem Confounding**: Verwendet genetische
  Instrumente, nicht nur beobachtete Exposition.
- **Respektiert die klinische Hierarchie**: Tod wird gegenüber
  Krankenhauseinweisung priorisiert, die gegenüber
  Biomarker-Verschlechterung priorisiert wird.
- **Berücksichtigt konkurrierende Risiken auf natürliche Weise**:
  Differentielles Überleben wird in den Estimanden aufgenommen und nicht
  wegkonditioniert.
- **Einzelnes zusammenfassendes Maß**: Eine DS-CWR ersetzt mehrere
  komponentenweise Analysen.
- **Diagnostisches Rahmenwerk**: SDPD testet die Exklusionsrestriktion;
  Q-Statistik testet Heterogenität; strukturierte Warnungen markieren
  Probleme.
- **Reproduzierbar**: Der Bootstrap ist Seed-kontrolliert; alle Ausgaben sind
  bei gegebenem Seed deterministisch.

### Limitationen

- **Daten auf Individualebene erforderlich**: Das Paket benötigt derzeit
  Daten mit einer Zeile pro Teilnehmer mit Genotyp, Exposition und
  Ergebnissen. MR nur mit Summary-Daten wird nicht unterstützt.
- **Nicht-Kollabierbarkeit**: Die Win-Ratio ist nicht kollabierbar. Die
  marginale DS-CWR unterscheidet sich von bedingten Win-Ratios. Dies ist eine
  Eigenschaft des Estimanden, kein Fehler.
- **Lineares PRS-Instrument**: Das Paket setzt einen linearen polygenetischen
  Risikoscore voraus. Nichtlineare Gen-Gen-Interaktionen werden nicht
  modelliert.
- **Diagonale GWAS-Kovarianz**: Der aktuelle Bootstrap verwendet
  SNP-für-SNP-Standardfehler, nicht die vollständige
  LD-Kovarianzmatrix. Dies kann die Unsicherheit unterschätzen, wenn
  SNPs im Kopplungsungleichgewicht stehen.
- **Pleiotropie-Anfälligkeit**: Die cCWR ist empfindlich gegenüber
  hierarchiekontaminierender Pleiotropie. Mortalitätsebenen-Pleiotropie
  von nur gamma = 0,05 kann die Abdeckung auf 12 % reduzieren.
- **Stichprobengrößenanforderungen**: Biobank-große Stichproben
  (N > 100.000) sind eine strenge statistische Voraussetzung für
  zuverlässige Inferenz.
- **Stratenzahl `D`**: Der Dezil-artige Schätzer erfordert die Wahl einer
  Stratenzahl. Ein kontinuierlicher, bandbreitengesteuerter Gradientenschätzer
  (von dem der Dezil-Schätzer der exakte Boxcar-Spezialfall ist) wurde
  implementiert und validiert, doch externe Prüfungen zeigten, dass er
  verrauschter als das Dezil ist und die `D`-Empfindlichkeit nicht reduziert
  (der instrumentenstandardisierte Gradient ist ein Quotient, und feinere
  Glättung lässt seinen Nenner schrumpfen). Die empfohlene Praxis ist der
  diskrete Schätzer mit einer Sensitivitätsanalyse über `D`, nicht eine
  kontinuierliche Reparametrisierung (siehe `dev/wp15-continuous-isg.md`).

### Leistungs- und Inferenz-Backends (optional)

Das Standard-Backend ist unverändert, aber die Arbeit der Phase II
(`dev/acceleration-roadmap.md`) fügte validierte, optionale Alternativen
über `mrwin_controls()` hinzu:

- `backend = "fast"` — ein subquadratischer, kompilierter (Rcpp) Win/Loss-Kernel
  (`Theta(N log^{K-1} N)` für `K` Prioritätsebenen), der bit-für-bit identisch
  mit dem dichten Backend ist und den Schätzer auf Biobank-Größe bringt (z. B.
  `K = 3`, `N = 80.000` in unter einer Sekunde).
- `inference = "analytic"` — eine geschlossene Influenzfunktions-Varianz, die den
  Multiplikator-Bootstrap reproduziert (mit einem exakten Monte-Carlo-Term für
  die Unsicherheit der GWAS-Gewichte) und die Bootstrap-Schleife für die
  Stichprobenkomponente entfernt.

> **Doppelt geordnete Stratifizierung — nicht empfohlen (negativer Befund).**
> `mrwin_doubly_ranked_strata()` (Tian/Burgess) ist implementiert, aber eine
> externe Kalibrierung zeigte, dass `stratification = "doubly_ranked"` mit dem
> Zwischen-Strata-DS-CWR-Estimand **inkompatibel** ist: Es balanciert das
> Instrument über die Strata hinweg, sodass der Kontrast benachbarter Strata
> konfounder-getrieben wird (Typ-I-Fehler ~1,0 unter Konfundierung). `mrwin()`
> führt es zwar weiterhin aus, gibt aber eine `doubly_ranked_invalid`-Warnung
> aus; der Standard ist `"prs_rank"`. Siehe `inst/spec/validation-findings.md`.

#### Konfidenzintervalle

Das **primär berichtete Intervall ist das Fieller-Intervall**. Eine externe
Kalibrierungssimulation zeigte, dass das ursprüngliche
Delta-Methoden-Quotientenintervall zu weit ist (es erreicht nicht den nominalen
Typ-I-Fehler), während die Fieller-Konstruktion auf derselben Kovarianz korrekt
kalibriert ist; das Delta-Methoden-Intervall wird nur als beschriftete Referenz
beibehalten. Bei einem schwachen Instrument wird das Fieller-Intervall ehrlich
als unbeschränkt berichtet, anstatt als fälschlicherweise enge Menge.

---

## Annahmen

Das `mrwin`-Rahmenwerk stützt sich auf die folgenden Annahmen:

1. **Relevanz**: Der PRS muss mit der Exposition assoziiert sein. Schwache
   Instrumente erzeugen instabile Schätzungen (prüfen Sie die
   `weak_instrument`-Warnung).

2. **Unabhängigkeit**: Der PRS darf nicht mit Confoundern assoziiert sein.
   Dies wird durch Mendelsche Segregation in homogenen Populationen
   erfüllt. Populationsstratifizierung kann dies verletzen; verwenden Sie
   Hauptkomponenten als Kovariaten.

3. **Exklusionsrestriktion**: Der PRS darf das Ergebnis nur über die
   Exposition beeinflussen. Die SDPD testet dies für den Endpunkt mit der
   höchsten Priorität. Bei Ablehnung kann die Schätzung ungültig sein.

4. **Monotonie**: Das Instrument muss die Exposition für alle Personen in
   dieselbe Richtung verschieben. Verletzungen können die Schätzung
   verzerren.

5. **Keine Effektmodifikation durch Überleben**: Der kausale Effekt sollte
   nicht zwischen Überlebenden und Nicht-Überlebenden variieren. Dies ist
   nicht testbar.

6. **Korrekte Prioritätsordnung**: Die klinische Hierarchie muss vor der
   Analyse festgelegt werden. Das Paket leitet sie nicht aus den Daten ab.

---

## Warncodes

| Code | Bedeutung | Was zu tun ist |
|---|---|---|
| `weak_instrument` | Phänotypische Verschiebung nahe Null oder Fieller-KI unbeschränkt | Schätzung als instabil betrachten |
| `sdpd_rejected` | MR-Egger-Achsenabschnitt deutet auf Pleiotropie hin | Sensitivitätsdiskussion hinzufügen |
| `sdpd_underpowered` | Zu wenige SNPs für SDPD | Power-Limitation berichten |
| `positivity_failure` | IPTW-ESS in einem Stratum zu klein | Kovariaten-Überlappung prüfen |
| `bridged_strata` | Strata aufgrund von Positivitätsversagen übersprungen | Überbrückungsschema berichten |
| `discordant_components` | Komponenten-MR-Richtungen widersprechen sich | Globalen Effekt nicht überbewerten |

---

## Simulations-Engine

Führen Sie Szenariogitter aus, um die Methode zu validieren:

```r
scenarios <- mrwin_scenarios(
  mrwin_config(n_outcome = 200, m_snps = 10, seed = 10),
  scenarios = c("A_null", "B_valid_IV", "C_pleiotropy", "D_hierarchy_discordant")
)

grid <- mrwin_run_simulation_grid(
  scenarios = scenarios, n_iter = 2,
  controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 11)
)

mrwin_simulation_summary(grid)
```

| Szenario | Bedeutung | Erwartetes Verhalten |
|---|---|---|
| A_null | Kein kausaler Effekt | DS-CWR nahe 1, Ablehnungsrate nahe 5 % |
| B_valid_IV | Valider kausaler Effekt | DS-CWR in erwarteter Richtung |
| C_pleiotropy | Direkte Pleiotropie | SDPD-Ablehnung nimmt zu |
| D_hierarchy_discordant | Entgegengesetzte Komponenteneffekte | Warnung wird ausgegeben |

---

## Dieses Paket zitieren

Wenn Sie `mrwin` in Ihrer Forschung verwenden, zitieren Sie bitte:

> Nguyen Thien Minh, N. Ahmad Aziz. „Causal Win Statistics: Integrating
> Instrumental Variable Estimation with Hierarchical Composite Endpoints."
> *arXiv preprint*, 2026.

BibTeX:

```bibtex
@article{NguyenAziz2026,
  title={Causal Win Statistics: Integrating Instrumental Variable
         Estimation with Hierarchical Composite Endpoints},
  author={Nguyen Thien Minh and N. Ahmad Aziz},
  journal={arXiv preprint},
  year={2026}
}
```

**Autoren:**
- Nguyen Thien Minh, Universität für Medizin und Pharmazie in Ho-Chi-Minh-Stadt,
  Vietnam (minhnt@ump.edu.vn)
- N. Ahmad Aziz (korrespondierender Autor), Deutsches Zentrum für
  Neurodegenerative Erkrankungen (DZNE) und Universität Bonn, Deutschland
  (Ahmad.Aziz@dzne.de)

---

## Verwandte Referenzen

**Win-Statistiken:**
- Pocock SJ, et al. The win ratio: a new approach to the analysis of
  composite endpoints in clinical trials. *Eur Heart J*. 2012;33(14):1744-1749.
- Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a
  composite endpoint based on prioritized components. *Biostatistics*.
  2016;17(1):178-191.
- Even Z, Josse A. Causal win ratio. *arXiv preprint*. 2025.

**Mendelsche Randomisierung:**
- Lawlor DA, et al. Mendelian randomization: using genes as instruments for
  making causal inferences in epidemiology. *Stat Med*. 2008;27(8):1133-1163.
- Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for
  causal inference in epidemiological studies. *Hum Mol Genet*.
  2014;23(R1):R89-R98.
- Bowden J, et al. Mendelian randomization with invalid instruments: effect
  estimation and bias detection through Egger regression. *Int J Epidemiol*.
  2015;44(2):512-525.

**Polygenetische Risikoscores:**
- Choi SW, et al. Tutorial: a guide to performing polygenic risk score
  analyses. *Nat Protoc*. 2020;15(9):2759-2772.

**Multi-State-Modelle:**
- Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks
  and multi-state models. *Stat Med*. 2007;26(11):2389-2430.

---

## Paketstruktur

```
mrwin/
  R/
    api.R              # High-Level-mrwin()-Workflow
    kernel.R           # Hierarchischer Vergleichskernel
    estimate.R         # DS-CWR-Schätzer, GLS-Pooling
    bootstrap.R        # Multiplikator-Bootstrap-Inferenz
    adjustment.R       # Ordinal-IPTW, ESS, Überbrückung
    sdpd.R             # Step-Down-Pleiotropie-Diagnostik
    simulate.R         # Simulations-Datengenerierungsprozess
    simulation_engine.R # Szenariogitter-Runner
    validation.R       # Eingabevalidierung
    methods.R          # print, summary, plot, tidy, report
    benchmark.R        # Leistungsbenchmarking
    backend_sparse.R   # Sparse-Backend
    config.R           # Simulationskonfiguration
    strata.R           # PRS-Stratumzuweisung
  tests/
    testthat/          # 273 Unit-Tests
  vignettes/           # 3 Vignetten
  inst/spec/           # Implementierungsspezifikationen
```

---

## Lizenz

MIT
