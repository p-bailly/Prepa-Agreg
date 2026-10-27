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
- Carte neuve proposée seulement si ses prérequis sont au moins « en cours » (ou sous le niveau de départ) ; une carte
  introduite dans la même séance compte. Ordre : niveau le plus bas, puis thèmes alternés, puis ordre du fichier.
- Un volet jamais noté d'une carte introduite (démonstration non ouverte) devient dû le lendemain de l'introduction.
- Séance : révisions d'abord (les plus en retard en premier), complétées par des neuves réparties entre elles ; au plus
  `seanceCartes` cartes et `seanceMinutes` minutes (la carte en cours est finie). Plafond de neuves par jour civil.
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
  "titre", "prevues" } ] } ] }`. `prevues` = nombre de cartes visées pour le chapitre (le total que vise le tableau de bord).
  L'ordre des chapitres d'un thème est l'ordre du parcours ; les niveaux y sont croissants.
- `content/cartes/<theme>/<niveau>.json` : `{ "cartes": [ … ] }`, l'ordre du fichier est l'ordre d'introduction dans un niveau.
- `content/exercices/<theme>/<niveau>.json` : `{ "exercices": [ … ] }`.
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
- Au moins 2 exercices de difficultés distinctes par carte, jusqu'à 6 pour les notions piégeuses ; jamais d'exercice de remplissage.
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
