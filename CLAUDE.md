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
