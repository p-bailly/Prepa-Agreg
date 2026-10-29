/// <reference types="vitest/config" />
import { readFileSync } from 'node:fs'
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// Numéro de version, le même que celui de l'APK (android/app/build.gradle) : « majeur.mineur » de package.json, puis le
// numéro de build de la CI (0 en local).
const [majeur, mineur] = (JSON.parse(readFileSync('package.json', 'utf8')) as { version: string }).version.split('.')
const version = `${majeur}.${mineur}.${process.env.GITHUB_RUN_NUMBER ?? 0}`

// Chemins relatifs : l'appli est servie depuis les fichiers de l'APK (Capacitor), sans serveur, et la même construction
// marche dans le sous-dossier de GitHub Pages (/Prepa-Agreg/).
export default defineConfig({
  base: './',
  plugins: [react()],
  // Numéro de version et date de construction, affichés dans les réglages pour savoir quelle version tourne.
  define: { __VERSION__: JSON.stringify(version), __BUILD__: JSON.stringify(new Date().toISOString().slice(0, 10)) },
  build: { chunkSizeWarningLimit: 1500 },
  test: {
    include: ['src/**/*.test.ts'],
    // Pas de cache des modules transformés : le résultat d'import.meta.glob (chargement de content/) y serait figé et un
    // fichier de contenu ajouté resterait invisible aux tests.
    fsModuleCache: false,
    testTimeout: 30000,
  },
})
