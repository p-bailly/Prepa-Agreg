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
