# CLAUDE.md

Prépa Agreg : appli Android perso (Capacitor + React 19 + Vite + TypeScript + KaTeX) pour monter de la Terminale jusqu'aux
attendus de l'écrit de l'agrégation externe de mathématiques : cartes de révision à répétition espacée, exercices de difficulté
croissante par notion, tableau de bord. Pas de version web publiée, pas de compte, pas de serveur. Code, commentaires, textes,
commits : en français. Les commentaires qui renvoient à une rubrique de ce fichier citent son titre : le garder en cas de renommage.

## Commandes
- `npm run dev` (5173, pour développer dans un navigateur ; l'appli n'est pas publiée en web), `npm test` (Vitest),
  `npm run typecheck`, `npm run build` (tsc + vite). Node 22.
- `python3 verification/tout.py` : vérifications calculatoires du contenu (SymPy, hors APK ; voir « Contrôle de l'exactitude »).
- APK en local : Android SDK (plateforme 36) + Java 21, `npm run android` → `android/app/build/outputs/apk/release/app-release.apk`.
  Sans Android Studio : `sdkmanager "platforms;android-36" "build-tools;36.0.0"` puis `sdk.dir=…` dans `android/local.properties`.
- CI : `verification.yml` (PR : typecheck, tests, SymPy, build, APK + contrôle du `versionName`) et `android.yml` (push sur la
  branche par défaut : tests, SymPy, build, APK publié dans la release « apk »).
- Logo : `npm run logo` (= `node design/logo/generer.mjs` puis `node design/logo/png.cjs`, qui pilote Chrome ; chemin dans
  `CHROME_PATH` au besoin). Les fichiers produits (SVG, calques adaptatifs, PNG) ne s'éditent jamais à la main.

## Android, versions, livraison (ne jamais changer : « Mes applis » en dépend)
- `applicationId` **`fr.perso.prepaagreg`**, APK **`prepa-agreg.apk`**, dépôt public `p-bailly/Prepa-Agreg`, nom « Prépa Agreg ».
  Lien permanent : `https://github.com/p-bailly/Prepa-Agreg/releases/latest/download/prepa-agreg.apk`. « Mes applis »
  (`p-bailly/Managing-Apps`) suivra `TrackedApp("Prépa Agreg", "p-bailly/Prepa-Agreg", "fr.perso.prepaagreg", "prepa-agreg.apk")`
  (ajout à faire par le propriétaire, pas ici).
- Signature : **la même clé que Scrabble et les autres applis**, copiée telle quelle de `p-bailly/Scrabble`
  (`android/app/scrabble.keystore`, alias et mots de passe dans `android/app/build.gradle`). Clé versionnée exprès, usage perso.
  Ne jamais en créer une autre : Android refuserait la mise à jour. Empreinte SHA-256 du certificat : `0B:D9:37:5B:…:85:17:8C`.
- Version `majeur.mineur.build` : `majeur.mineur` = champ `version` de `package.json` (à relever à la main), `build` =
  `GITHUB_RUN_NUMBER` (0 en local) ; `versionCode` = `github.run_number`. Utilisée par l'APK, `vite.config.ts` (`__VERSION__`)
  et les Réglages ; `verification.yml` contrôle l'APK.
- Polices Latin Modern : `.otf` CTAN non modifiés dans `src/polices/` (GUST Font License, voir `src/polices/LISEZMOI.md`).
  Police d'interface Plus Jakarta Sans : paquet npm `@fontsource-variable/plus-jakarta-sans` (SIL OFL 1.1), importé dans
  `src/main.tsx`, donc embarqué par Vite.

## Architecture
- Logique pure en TypeScript, sans React, testée avec Vitest et **une horloge injectée** (`src/lib/temps.ts` : `Horloge`,
  `jourLocal(t, décalage)`) ; jamais de `Date.now()` hors de `horlogeSysteme`.
  - `src/repetition/` : `planification.ts` (SM-2 simplifié par volet), `etat.ts` (neuve / en cours / acquise, prérequis),
    `seance.ts` (cartes à revoir, cartes neuves, plafond quotidien, filtres, ordre d'introduction), `modes.ts` (quatre modes),
    `deroulement.ts` (séance courte : reprises, arrêt au temps).
  - `src/exercices/niveaux.ts` : cases des cinq niveaux, exercice suggéré. `src/stats/tableau.ts` : tableau de bord.
  - `src/contenu/` : types, chargement (`import.meta.glob` de `content/`), découpage du texte (`texte.ts`), validation.
  - `src/lib/` : stockage local, réglages, progression (fonctions pures), sauvegarde JSON, notifications, partage de fichier.
- Écrans React dans `src/pages/` (Accueil, Parcours, Chapitre, Stats, Réglages, Séance, EcranCarte, EcranExercice) ; données
  partagées par `src/etat/Donnees.tsx` (contexte : contenu, réglages, progression, horloge) ; navigation dans `App.tsx`
  (quatre onglets + pile d'écrans, bouton retour Android). Dossiers prévus pour plus tard (pas de code en v1) :
  `src/sujets/`, `src/seance/` (sujets chronométrés, séances longues) et `content/sujets/`.
- Persistance locale versionnée : `reglages-v1`, `progression-v1` (localStorage de la WebView) ; changer un format = changer le
  suffixe ; la lecture fusionne avec les valeurs par défaut et répare les champs invalides (`normaliser…`). La progression est
  indexée par identifiant de carte (`idCarte#enonce`, `#demo`, `#contre`) et d'exercice : **renommer un identifiant perd son
  historique**, ne jamais le faire sans migration.
- Tout est embarqué dans l'APK (contenu, KaTeX et ses polices, Latin Modern) : hors ligne par construction, aucun CDN.

## Allure visuelle
Thème « Tonal ardoise » (maquette 7b), tout dans `src/styles.css` sous forme de jetons (`--fond`, `--surface`, `--primaire`,
`--conteneur`, `--sable`, `--menthe`, `--rose`, `--lilas`…), redéfinis pour le sombre sous `:root[data-theme='dark']`.
- Fond gris bleuté `#F4F6FA`, surfaces blanches à grands arrondis (20 à 32 px), bleu ardoise `#3E5F8A`, aplats tonals ; ni
  ombre ni dégradé. Sombre : fond `#111418`, mêmes rôles (`index.html` et `capacitor.config.json` reprennent les fonds).
- Interface en Plus Jakarta Sans ; textes du contenu (`.texte` : énoncés, démonstrations, corrections) en Latin Modern, pour
  aller avec les formules KaTeX.
- Accueil : bloc « Aujourd’hui » bleu qui contient le bouton Réviser, tuiles sable / menthe / rose, objectif du jour.
  Carte : bandeau tonal (thème, niveau, titre), encadrés en aplat selon le type (théorème bleu, définition menthe, exemple
  lilas, méthode sable, remarque neutre). Notation : quatre aplats collés en pilule (rose, sable, bleu, menthe), chacun avec
  son nom et son délai (jamais la couleur seule). Onglets : icône dans une pilule pour l'onglet actif.
- Le logo (crayon rouge `#CE2029`) n'a pas suivi : il se régénère avec `npm run logo` si on veut l'accorder.

## Répétition espacée et séance (choix non évidents)
- Notes : À revoir (reprise dans 10 min, la carte revient en fin de séance, deux fois au plus), Difficile, Bien, Facile.
  Première fois : 1 j, 4 j, 9 j ; après un oubli : 1 j, 2 j, 4 j ; ensuite Difficile = max(I+1, 1,2 I), Bien = (I + retard/2)
  × facilité, Facile = (I + retard) × facilité × 1,3, strictement croissants, plafond 365 j ; facilité 1,3 à 3. Échéances en
  **jours civils** (un changement de fuseau ne décale rien de plus d'un jour). Volet **acquis** à 21 j ; carte acquise quand
  tous ses volets le sont.
- **On apprend avant de réviser.** Une notion nouvelle se découvre dans une **leçon** (`EcranLecon` : introduction du
  chapitre `content/lecons/<chapitre>.json`, puis par notion le texte `cours`, l'énoncé, l'exemple, la démonstration, les
  contre-exemples, les pièges, les exercices) ; passer à la suite la marque apprise (`progression.apprises`, carte → jour).
  Une carte n'entre en révision qu'apprise, et jamais le jour même : tous ses volets deviennent dus le lendemain.
- Leçon du jour : les prochaines notions à apprendre (prérequis au moins appris, ou sous le niveau de départ ; ordre :
  niveau le plus bas, thèmes alternés, chapitres de `themes.json`, ordre du fichier), dans la limite du plafond quotidien
  (`plafondNeuves`, compté par jour civil), regroupées par chapitre. Accueil : « Apprendre » est le bouton principal quand
  rien n'est à revoir. Parcours → chapitre → « Apprendre / Relire la leçon » à tout moment.
- Séance de révision : cartes apprises dues, les plus en retard en premier ; au plus `seanceCartes` cartes et
  `seanceMinutes` minutes (la carte en cours est finie).
- Modes : une carte neuve se découvre toujours en « énoncer de tête » ; sinon on prend le mode possible et actif qui suit le
  dernier utilisé (`progression.dernierMode`). Volet principal = premier volet dû (énoncé, démonstration, contre-exemple), noté
  par le pied ; les autres volets dus ont leurs boutons dans le volet ouvert.
- « Niveau atteint » : plus haut niveau dont 80 % des cartes **visées** (`prevues` de `themes.json`, au moins les cartes
  rédigées) sont acquises ; les niveaux sous le niveau de départ comptent comme acquis (ou sont masqués, au choix).
- Exercices : hors répétition espacée ; niveau ■ réussi quand chacun de ses exercices a été réussi une fois, ◧ dès une
  tentative ; exercice suggéré = premier non réussi du premier niveau non réussi.

## Niveaux (choix, pas une source officielle)
Plancher : programme de Terminale (spécialité mathématiques et option mathématiques expertes, programmes de 2019 en vigueur en
2026-2027 ; un nouveau programme de spécialité, BO n° 14 du 2 avril 2026, s'appliquera à la rentrée 2027). Plafond : programme
de l'écrit de l'agrégation externe (session 2027 : douze titres, de l'algèbre linéaire aux méthodes numériques ; épreuves
« composition de mathématiques générales » et « composition d'analyse et probabilités », 6 h, coefficient 4 chacune, d'après
devenirenseignant.gouv.fr consulté le 2026-10-27). **Entre les deux (L1, L2, L3), il n'y a pas de programme officiel** : le
découpage de `themes.json` est construit à rebours depuis le plafond et c'est un choix, pas une source officielle.

## Format du contenu

Le contenu est fait de **fichiers de données** : ajouter du contenu = ajouter un fichier, jamais toucher au code. Tous les
fichiers de `content/` sont chargés par `src/contenu/chargement.ts` (`import.meta.glob`, embarqués dans l'APK) et contrôlés par
`src/contenu/validation.test.ts`.

- `content/themes.json` : la carte du parcours complet. `{ "themes": [ { "id", "titre", "chapitres": [ { "id", "niveau",
  "titre", "notions", "prevues" } ] } ] }`. `notions` = les notions à couvrir, une carte chacune (plan de rédaction, établi
  d'après les sous-sections du programme de l'agrégation 2027 et des programmes de Terminale) ; `prevues` = leur nombre (le
  total que vise le tableau de bord ; contrôlé par le test). L'ordre des chapitres d'un thème est l'ordre du parcours, les
  niveaux y sont croissants.
- `content/cartes/<theme>/<niveau>.json` ou `<niveau>-<chapitre>.json` (fichier réservé à un chapitre, pour rédiger des
  chapitres en parallèle) : `{ "cartes": [ … ] }`. Ordre d'introduction : chapitres dans l'ordre de `themes.json`, puis ordre
  du fichier.
- `content/exercices/<theme>/<niveau>.json` ou `<niveau>-<chapitre>.json` : `{ "exercices": [ … ] }`.
- Niveaux, du plancher au plafond : `terminale`, `l1`, `l2`, `l3`, `agreg`.

### Carte
`{ "id", "theme", "niveau", "chapitre", "type", "titre", "recto", "hypotheses", "hypotheses_pieges", "enonce", "enonce_trous",
"trous", "demo", "exemple", "contre_exemple", "affirmations", "prerequis", "doute" }`
- Obligatoires : `id` (unique, en minuscules avec tirets, **jamais renommé** : la progression y est attachée), `theme`, `niveau`,
  `chapitre` (un `id` de chapitre du thème, de même niveau), `type` (`theoreme`, `definition` ou `methode`), `titre`, `recto` (la
  consigne, ex. « Énoncer le théorème de Rolle avec ses hypothèses. »), `enonce` (le résultat complet, hypothèses comprises),
  `prerequis` (identifiants de cartes, éventuellement `[]`).
- Facultatifs : `hypotheses` et `hypotheses_pieges` (listes de phrases ; le mode « hypothèses à cocher » demande les deux),
  `enonce_trous` + `trous`, `demo`, `exemple`, `contre_exemple`, `affirmations`, `doute`.
- `enonce_trous` : l'énoncé où des morceaux sont remplacés par `{{1}}`, `{{2}}`… (jamais à l'intérieur d'une formule `$…$`) ;
  `trous[i]` = `{ "reponse", "leurres": [3 textes] }` correspond à `{{i+1}}`. Les leurres sont crédibles et **faux** (l'énoncé
  obtenu avec un leurre est faux ou n'a pas de sens).
- `affirmations` : `[{ "texte", "vrai": true|false, "contre_exemple", "justification" }]` ; une affirmation fausse a un
  `contre_exemple` (vérifié) ; une vraie peut avoir une `justification`. Au moins une fausse si le champ est présent.
- Volets notés séparément en répétition espacée : `enonce` (toujours), `demo` et `contre` (si `demo` / `contre_exemple`).
- `doute` : présent seulement si une partie reste incertaine ; l'appli affiche « Non garanti : … ».

### Exercice
`{ "id", "carte", "difficulte", "titre", "enonce", "duree_min", "indices", "reponse_finale", "correction", "doute" }`
- `difficulte` de 1 à 5, relative à la notion : 1 application directe ; 2 calcul ou application standard ; 3 classique en plusieurs
  étapes ; 4 question d'écrit (enchaînement, hypothèses à vérifier) ; 5 difficile (idée non immédiate, plusieurs notions).
- Par carte : 3 à 6 exercices progressifs, au moins un de difficulté 1 ou 2 et au moins un de difficulté 4 ou 5 (contrôlé par
  le test) ; aux niveaux L3 et agrégation, un problème de difficulté 5 en questions enchaînées quand c'est pertinent. Jamais
  d'exercice de remplissage.
- `variantes` (facultatif) : autres jeux de valeurs, `[{ "titre", "enonce", "reponse_finale", "correction", "indices" }]`
  (`reponse_finale` obligatoire si l'exercice en a une ; `indices` obligatoires si l'exercice en a, contrôlé par le test).
  Une variante sans `titre` reprend celui de l'exercice, qui ne doit donc citer aucune valeur propre à l'exercice. Tout exercice à valeurs numériques en a deux. L'appli
  propose une variante jamais faite, sinon la plus anciennement faite ; la tentative garde son numéro de variante (0 = exercice
  tel qu'écrit). Ne jamais réordonner les variantes d'un exercice déjà publié (l'historique y renvoie par numéro).
- `indices` : 2 ou 3 à partir de la difficulté 3 (idée clé, puis première étape) ; ils ne donnent pas la réponse.
- `reponse_finale` (seulement pour un exercice à résultat : calcul, valeur, vrai/faux) : `{ "question", "propositions": [4 textes
  distincts], "bonne": indice 0..3 }` ; l'appli mélange les propositions.
- `correction` : rédigée comme à l'écrit.

### Texte et formules
- Texte avec formules KaTeX : `$…$` en ligne, `$$…$$` centré. Paragraphes séparés par une ligne vide (`\n\n`), lignes de liste
  commençant par `- `, gras `**…**`. Pas d'autre balisage. En JSON, chaque `\` s'écrit `\\`.
- Intervalles à la française : `\left]a,b\right[`, `\left[a,b\right[`. Ensembles `\mathbb{R}`, `\mathbb{N}`. Pas de `\dfrac` dans
  les formules en ligne trop hautes ; `\displaystyle` à éviter.
- Toute formule doit passer `katex.renderToString(…, { throwOnError: true })` (contrôlé par le test).

## Contrôle de l'exactitude (personne ne relit le contenu)
Un énoncé faux appris par cœur est pire qu'une carte absente : mieux vaut peu de contenu sûr que beaucoup de contenu douteux.
- Cartes : (1) rédaction ; (2) relecture par un **sous-agent indépendant** qui ne voit pas le raisonnement de l'auteur et
  re-dérive énoncé, hypothèses, pièges, trous (chaque leurre doit rendre l'énoncé faux), démonstration, contre-exemples et
  affirmations ; (3) correction ou retrait de tout ce qu'il conteste ; (4) faits calculatoires vérifiés dans
  `verification/<theme>.py` (SymPy, `assert` commentés citant l'identifiant ; `python3 verification/tout.py`, lancé par la CI).
- Exercices : un sous-agent indépendant **résout d'abord chaque exercice à partir des seuls énoncés** (sans corrections, indices
  ni propositions), écrit ses réponses, puis compare ; tout écart est tranché par un calcul (SymPy) ou une preuve détaillée. Il
  vérifie aussi énoncé complet et non ambigu, bonne proposition seule juste, leurres faux, indices qui ne donnent pas la réponse,
  difficulté plausible.
- Un élément dont une partie reste incertaine est retiré, ou gardé avec un champ `doute` (bandeau « Non garanti : … »).
- Chaque correction de contenu après un signalement (Réglages → Erreurs signalées, inclus dans l'export) suit la même procédure.

### Journal
- 2026-10-27 · `cartes/arithmetique/terminale.json`, `cartes/arithmetique/l1.json` (6 cartes) : relecture indépendante, tout
  re-dérivé, calculs en Python. Corrigé : exemple de division d'un négatif (cas `b | a`), piège de Fermat qui découlait des vraies
  hypothèses, piège faible du lemme chinois remplacé (« m divise n »), prérequis du lemme chinois (+ Gauss).
- 2026-10-27 · `exercices/arithmetique/terminale.json`, `exercices/arithmetique/l1.json` (22 exercices) : résolution à
  l'aveugle, toutes les réponses concordent, leurres vérifiés faux par force brute. Corrigé : trois difficultés surestimées
  (division ex2 3→2, ex3 4→3, Bézout ex3 4→3), question ambiguë sur la « réciproque de Fermat » (Fermat ex4), indice qui
  éliminait un leurre (Gauss ex4). `verification/arithmetique.py` (dont 2^340 ≡ 1 mod 341, X² − 1 dans Z/8Z).
- 2026-10-27 · `cartes/probabilites/terminale.json`, `cartes/probabilites/l1.json` (5 cartes) : relecture indépendante, calculs
  exacts en fractions. Corrigé : remarque de la démonstration de la loi faible (« non presque sûrement constantes »), sens de
  P(B | A_i) précisé (Bayes), n ≥ 3 dans le texte à trous de l'indépendance mutuelle, prérequis inutiles retirés.
- 2026-10-27 · `exercices/probabilites/terminale.json`, `exercices/probabilites/l1.json` (19 exercices) : résolution à l'aveugle,
  19 réponses concordantes, leurres vérifiés faux en calcul exact. Corrigé : notation « a = 2 » non introduite (Bienaymé ex3).
  `verification/probabilites.py` (dont les deux pièces A, B, C deux à deux indépendants, cas d'égalité de Bienaymé-Tchebychev).

- 2026-10-27 · `cartes/analyse-reelle/terminale.json`, `cartes/analyse-reelle/l1.json`, `exercices/analyse-reelle/terminale.json`, `exercices/analyse-reelle/l1.json`, `lecons/ar-t-{suites,continuite,derivation}.json`, `lecons/ar-l1-{reels,suites,continuite,derivabilite}.json` (10 cartes, 36 exercices, 42 variantes, 7 introductions) : relecture complète indépendante (cartes, `cours` et introductions re-dérivés ; exercices et variantes résolus à l'aveugle, 36 réponses concordantes, leurres vérifiés faux avec SymPy) ; corrigé : 4 pièges qui découlaient des vraies hypothèses (gendarmes ×2, TVI, bornes atteintes), prérequis du TVI (+ limite monotone), phrases de cours (point fixe avec f décroissante, minimum seulement pour bornes atteintes), introductions (continuité « sans lever le crayon » limitée à un intervalle, hypothèses du point fixe, définition de Cauchy), triple saut de ligne (limite-monotone-ex2), indice qui livrait le contre-exemple (bornes-atteintes-ex3) ; aucun doute.
- 2026-10-27 · content/exercices/arithmetique/terminale.json, l1.json, content/cartes/arithmetique/terminale.json, l1.json (champ cours), content/lecons/arith-t-divisibilite, arith-t-pgcd, arith-t-premiers, arith-l1-congruences, arith-l1-polynomes.json, verification/arithmetique.py (6 cartes, 24 exercices dont 2 ajoutés, 44 variantes) : la relecture a résolu à l'aveugle les 2 exercices ajoutés et les 44 variantes (réponses toutes concordantes, leurres faux, ensembles de solutions comparés par force brute) et relu les 6 cours et les 5 introductions contre les programmes ; corrigé : 6 titres rendus génériques et indices propres aux variantes de division-ex1, phrase sur l'identité en Q4 de division-ex4, énoncé de fermat-ex4 v2 (étape modulo 19), hypothèses manquantes dans les cours (test de Fermat avec n ∤ a, points distincts, conclusion de Bézout) et les introductions (Euler avec a ∧ n = 1, b ≥ 1, points distincts, Wilson précisé) ; aucun doute.
- 2026-10-27 · `exercices/probabilites/terminale.json`, `exercices/probabilites/l1.json` (enrichissement : 1 exercice ajouté, proba-t-binomiale-ex4, et 36 variantes ; au total 20 exercices, 36 variantes), `cartes/probabilites/terminale.json`, `cartes/probabilites/l1.json` (5 cartes : champ `cours`), `lecons/proba-t-binomiale.json`, `lecons/proba-t-concentration.json`, `lecons/proba-l1-conditionnement.json` : exercice et variantes résolus à l'aveugle (37 réponses concordantes, leurres vérifiés faux en calcul exact), cours et introductions re-dérivés et confrontés au programme 2027 ; corrigé : intuition et piège « sans remise » du cours binomial (loi non B(n,p) dès n ≥ 2 et urne mixte), piège d'indépendance et notations du cours de la loi faible, « additivité de la variance » nuancée (covariances nulles) dans l'introduction, deux retouches de style, indices propres à binomiale-ex3 variante 1, indice adouci (indépendance-mutuelle-ex3 variante 2) ; `verification/probabilites.py` complété ; aucun doute.

- 2026-10-27 · `cartes/integration/terminale-int-t-integrale.json`, `exercices/integration/terminale-int-t-integrale.json`, `lecons/int-t-integrale.json`, `verification/int_t_integrale.py` (5 cartes, 23 exercices, 46 variantes) : relecture indépendante des cartes (tout re-dérivé en SymPy, 768 formules KaTeX, leurres et contre-exemples vérifiés) et des exercices (69 versions résolues à l'aveugle, 69 réponses concordantes, leurres vérifiés faux) ; corrigé : phrase du cours de la valeur moyenne (« sauf pour une fonction affine » → « en général »), indices propres à chacune des 46 variantes (dont 8 indices hérités inadaptés), titres propres pour 4 variantes, aire-ex5 en difficulté 4 au lieu de 5 ; doutes : les démonstrations admettent les propriétés de l'aire et le théorème des bornes atteintes (L1), et quelques difficultés 4 pourraient valoir 3 selon le relecteur.
- 2026-10-27 · `cartes/fondements/l1-fond-l1-ensembles.json`, `exercices/fondements/l1-fond-l1-ensembles.json`, `lecons/fond-l1-ensembles.json` (6 cartes, 29 exercices, 34 variantes) : relecture indépendante des cartes (énoncés, pièges, leurres, démonstrations, contre-exemples re-dérivés, faits finis vérifiés en Python) et résolution à l'aveugle des exercices (toutes les réponses concordent, 51 questions dont chaque leurre testé faux, formules KaTeX contrôlées), `verification/fond_l1_ensembles.py` ; corrigé : cours des tiroirs (finitude des tiroirs, pas des ensembles), des dénombrables (produits et réunions « peuvent ne pas être » au plus dénombrables) et de Cantor (sens de l'injection, chiffres 4/5 suffisants), texte à trous de Cantor (« par l'argument diagonal »), contre-exemples des dénombrables rendus indépendants de la carte Cantor, chiffres d5, d6… fixés (cantor-ex1), indices propres aux variantes de Zeckendorf, titres de variantes (timbres, dates de naissance), deux difficultés 5 → 4 (tiroirs-ex5, denombrables-ex5) ; aucun doute.
- 2026-10-27 · `cartes/equations-differentielles/terminale-ed-t-premier-ordre.json`, `exercices/equations-differentielles/terminale-ed-t-premier-ordre.json`, `lecons/ed-t-premier-ordre.json` (3 cartes, 15 exercices, 26 variantes) : relecture indépendante des cartes (énoncés, hypothèses, pièges, trous, démonstrations, contre-exemples et affirmations re-dérivés, SymPy, 569 formules KaTeX) et résolution à l'aveugle des 41 énoncés (réponses concordantes, chaque proposition testée par SymPy) ; corrigé : accord du trou 3 de `ed-t-premier-ordre-ay`, hypothèses de la forme u'u^n dans le cours des primitives, phrase fausse de l'introduction (vitesse dépendant « de façon affine » de la valeur), indices propres à toutes les variantes (dont trois indices hérités faux), énoncé à trous de `-ay-b` qui donnait les réponses, difficulté de `ay-ex5` 5→4, « l'ensemble des fonctions » dans `primitives-ex4` ; `verification/ed_t_premier_ordre.py` ; aucun doute.
- 2026-10-27 · `cartes/analyse-reelle/terminale-ar-t-suites.json`, `cartes/analyse-reelle/terminale-ar-t-continuite.json`, `exercices/analyse-reelle/terminale-ar-t-suites.json`, `exercices/analyse-reelle/terminale-ar-t-continuite.json`, `verification/ar_t_suites.py`, `verification/ar_t_continuite.py` (6 cartes, 26 exercices, 42 variantes) : relecture indépendante des cartes (énoncés, pièges, leurres, démonstrations, contre-exemples, affirmations et `cours` re-dérivés, 939 formules KaTeX, faits vérifiés avec SymPy) et des exercices (68 éléments résolus à l'aveugle, toutes les réponses concordantes, leurres vérifiés faux) ; corrigé : 4 pièges qui découlaient des vraies hypothèses réécrits en remplacements explicites (récurrence ×2, bijection ×2), `cours` de la suite géométrique ($u_0\neq0$, signe de $u_0$, $q\neq1$), prérequis de la suite géométrique (+ opérations), indice de operations-ex4 et de ses 2 variantes ($n^2$ en facteur sous la racine), titre générique pour image-suite-ex2 (variante sinus) ; aucun doute.
- 2026-10-27 · `cartes/algebre-lineaire/terminale-al-t-matrices.json`, `exercices/algebre-lineaire/terminale-al-t-matrices.json`, `lecons/al-t-matrices.json` (5 cartes, 23 exercices, 44 variantes) : relecture indépendante des cartes (tout re-dérivé en SymPy, 905 formules KaTeX, théorème des chemins testé par énumération) et des exercices (résolution à l'aveugle, 67 réponses concordantes, leurres vérifiés faux par script) ; corrigé : pièges qui découlaient des vraies hypothèses (produit, suites, Markov) et hypothèse du carré d'une somme ajoutée, piège du cours des suites (convergence « quelle que soit U_0 »), déterminant explicité dans suites-ex2 et variantes, « sommets donnés » dans adjacence-ex3, difficultés produit-ex3 3→2, produit-ex4 4→3, adjacence-ex5 5→4 ; `verification/al_t_matrices.py` ; aucun doute.
- 2026-10-27 · `cartes/geometrie/terminale-geo-t-espace.json`, `exercices/geometrie/terminale-geo-t-espace.json`, `lecons/geo-t-espace.json`, `verification/geo_t_espace.py` (5 cartes, 20 exercices, 38 variantes) : relecture indépendante des cartes (énoncés, hypothèses, trous, démonstrations, contre-exemples, affirmations, `cours` et introduction re-dérivés ; 925 formules KaTeX ; calculs refaits en SymPy) et des exercices (58 versions résolues à l'aveugle, toutes concordantes, leurres vérifiés faux, propositions d'équations de plans testées par proportionnalité) ; corrigé : 2 pièges de coplanaires qui découlaient des vraies hypothèses (« seulement deux à deux non colinéaires (à la place de non coplanaires) », « u et v orthogonaux »), phrases de cours fausses (deux équations « bien choisies » pour la décomposition, multiplication de a, b, c seulement si d ≠ 0, sphère de centre M et triangle MHN non dégénéré dans projete), précision « on fixe une autre inconnue » (equation-plan), prérequis droite-param ajouté à projete, « parallèles et distinctes » (coplanaires-ex4), non-colinéarités justifiées (orthogonalite-ex4 v2, projete-ex4 v1/v2), notation D_m (droite-param-ex4) ; doutes : la base de l'espace admet la propriété d'incidence droite/plan, la formule de distance point-plan dépasse un peu le programme, plusieurs difficultés pourraient valoir un cran de moins selon le relecteur.
- 2026-10-27 · `cartes/analyse-reelle/terminale-ar-t-derivation.json`, `exercices/analyse-reelle/terminale-ar-t-derivation.json`, `verification/ar_t_derivation.py`, retouche d'une phrase de `lecons/ar-t-derivation.json` (4 cartes, 18 exercices, 32 variantes) : relecture indépendante des cartes (énoncés, démonstrations, pièges, leurres et contre-exemples re-dérivés, dont $x^3+x^4\sin\frac1x$, calculs SymPy) et des exercices résolus à l'aveugle (50 réponses concordantes, leurres vérifiés faux) ; corrigé : trois phrases de cours (composée : « car » abusif, voir $\sqrt{x^4}$ ; inflexion : méthode qui excluait à tort les zéros de $f''$ sans changement de signe ; tangentes : $g'$ de signe large), introduction alignée sur la définition « la courbe traverse sa tangente », composee-ex4 (question 3 ouverte, le QCM n'est plus dans l'énoncé), convexite-ex4 moins guidé (difficulté 4 gardée), indice de inflexion-ex4 et variantes qui livrait le piège, tangentes-ex5 5 → 4 ; aucun doute (définition du point d'inflexion explicitée, l'autre convention signalée).
- 2026-10-27 · `cartes/probabilites/terminale-proba-t-binomiale.json`, `cartes/probabilites/terminale-proba-t-concentration.json`, `exercices/probabilites/terminale-proba-t-binomiale.json`, `exercices/probabilites/terminale-proba-t-concentration.json`, `verification/proba_t_binomiale.py`, `verification/proba_t_concentration.py` (5 cartes, 20 exercices, 36 variantes) : relecture indépendante des cartes (énoncés, hypothèses, pièges, leurres des trous, démonstrations, contre-exemples et affirmations re-dérivés, calculs exacts et force brute) et résolution à l'aveugle des exercices (56 réponses sur 56 concordantes, leurres vérifiés faux, propositions distinctes) ; corrigé : piège du cours de dénombrement (facteur k! seulement pour des éléments distincts), contre-exemple « n!/k! » de Pascal (5!/2! au lieu d'un cas où k = n − k), « une anagramme » (denombrement-ex3 et variantes), phrase de produit-ex2 (probabilités « pas toutes égales », fausse pour la variante 2) et titre neutre « Trois machines indépendantes », notation des indicatrices expliquée et « pour n ≥ 2 » dans la carte de linéarité ; aucun doute (les difficultés de produit-ex3 et denombrement-ex3 pourraient valoir un cran de moins selon le relecteur).
- 2026-10-27 · `cartes/geometrie/terminale-geo-t-complexes.json`, `exercices/geometrie/terminale-geo-t-complexes.json`, `lecons/geo-t-complexes.json`, `verification/geo_t_complexes.py` (6 cartes, 25 exercices, 48 variantes) : relecture indépendante des cartes (énoncés, hypothèses, pièges, trous, démonstrations, contre-exemples et affirmations re-dérivés en SymPy, 1024 formules KaTeX, introduction confrontée aux programmes) et résolution à l'aveugle des 73 énoncés (réponses toutes concordantes, leurres vérifiés faux en SymPy) ; corrigé : cours des racines de l'unité (polygone seulement pour n ≥ 3, « premier exemple » → « exemple frappant », z^n = a avec a ≠ 0), cours du second degré (rôle du signe de Δ) et démonstration (deux solutions réelles distinctes si Δ > 0), question de second-degre-ex3 et de ses variantes (ensemble des solutions de P(z) = 0), indices propres aux variantes d'exponentielle-ex2, indice de geometrie-ex4 qui livrait la réponse, et deux précisions facultatives (pourquoi le conjugué de w est non nul ; le point A appartient à F dans conjugue-ex4) ; aucun doute.
- 2026-10-27 · content/cartes/analyse-reelle/terminale-ar-t-fonctions-usuelles.json, content/exercices/analyse-reelle/terminale-ar-t-fonctions-usuelles.json, content/lecons/ar-t-fonctions-usuelles.json, verification/ar_t_fonctions_usuelles.py (5 cartes, 23 exercices, 40 variantes) : relecture indépendante des cartes (énoncés, pièges, leurres, démonstrations, 38 contrôles SymPy, 848 formules KaTeX) et des exercices (63/63 réponses retrouvées à l'aveugle, leurres et valeurs citées recalculés) ; corrigé : pièges de ln, croissances et limites-zero réécrits en remplacements explicites (ils découlaient des vraies hypothèses), limites de exp justifiées dans la démo de ln, continuité de cos en 0 justifiée dans la preuve de sin x/x, phrase fausse du cours des croissances sur les produits, contre-exemple « limite finie » remplacé par |x|, indice de trigo-ex3 V1 qui éliminait un leurre, signe de f' en 0 (croissances-ex3), coquille ln-ex2 V1 ; restent admises sans preuve : l'existence de solutions non continues de g(x+y)=g(x)g(y) (axiome du choix) et les formules d'addition de sin et cos.

## Pièges
- Formules : `\\` en JSON ; jamais de trou `{{n}}` dans une formule (le test le refuse) ; intervalles `\left]…\right[`.
- `content/` est importé à la construction : après un ajout de fichier, relancer `npm run dev`.
- Les symboles ■ ◧ □ ● ◐ ○, le cadenas, les icônes d'onglets et le chevron des volets sont dessinés en SVG
  (`src/composants/Icones.tsx`) : Latin Modern et Plus Jakarta Sans n'ont pas ces glyphes.
- Largeurs des barres : pourcentages CSS sans espace (`50%`), l'affichage français met une espace (`50 %`).
- Rappel quotidien : `@capacitor/local-notifications` (`schedule.on`), permission Android 13+ demandée à l'activation ;
  dans un navigateur, il n'y a pas de rappel.
- Export : fichier écrit dans le cache puis feuille de partage Android (`@capacitor/filesystem` + `@capacitor/share`) ; import :
  `<input type="file">` (sélecteur de fichiers Android). L'import remplace la progression après confirmation.
