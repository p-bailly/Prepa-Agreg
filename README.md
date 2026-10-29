# Prépa Agreg

Appli Android perso pour **monter de la Terminale jusqu'aux attendus de l'écrit de l'agrégation externe de
mathématiques** : cartes de révision à répétition espacée (énoncé, démonstration et contre-exemple notés séparément),
exercices de difficulté 1 à 5 pour chaque notion, tableau de bord de progression. Tout est sur le téléphone : pas de
compte, pas de serveur, hors ligne.

Contenu rédigé et contrôlé par une IA, **non relu par un humain** (méthode de contrôle dans [`CLAUDE.md`](CLAUDE.md)).

## Installer

Version web : <https://p-bailly.github.io/Prepa-Agreg/> (republiée à chaque push sur la branche par défaut ; la progression
reste dans le navigateur, export / import JSON dans les Réglages pour la transférer).

APK : <https://github.com/p-bailly/Prepa-Agreg/releases/latest/download/prepa-agreg.apk>. Reconstruit par GitHub Actions
à chaque push sur la branche par défaut, il s'installe par-dessus l'ancien et est suivi par
[« Mes applis »](https://github.com/p-bailly/Managing-Apps).

## Commandes

Prérequis : Node.js 22 (et Python 3 + SymPy pour les vérifications du contenu).

```bash
npm install
npm run dev               # http://localhost:5173, pour développer dans un navigateur
npm test                  # tests Vitest (logique, contenu, formules KaTeX)
npm run typecheck
npm run build
python3 verification/tout.py   # vérifications calculatoires du contenu (SymPy)
npm run android           # APK en local (Android SDK 36 + Java 21)
```

Architecture, formats du contenu, journal du contrôle de l'exactitude : [`CLAUDE.md`](CLAUDE.md).
