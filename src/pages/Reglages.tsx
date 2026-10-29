// Réglages : niveau de départ, niveaux actifs, thème, export et import de la progression.
import { useRef, useState } from 'react'
import { NIVEAUX, NOMS_NIVEAUX, type Niveau } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { partagerFichier } from '../lib/fichiers'
import type { Reglages as TypeReglages } from '../lib/reglages'
import { fabriquerSauvegarde, lireSauvegarde } from '../lib/sauvegarde'

export function Reglages() {
  const { reglages: r, changerReglages, progression, remplacer, horloge, aujourdhui, contenu } = useDonnees()
  const [message, setMessage] = useState<string | null>(null)
  const fichier = useRef<HTMLInputElement>(null)
  const maj = (f: Partial<TypeReglages>) => changerReglages((x) => ({ ...x, ...f }))
  const basculerListe = <T,>(l: T[], v: T): T[] => (l.includes(v) ? l.filter((x) => x !== v) : [...l, v])

  const exporter = async () => {
    try {
      const nom = `prepa-agreg-${aujourdhui()}.json`
      await partagerFichier(nom, JSON.stringify(fabriquerSauvegarde(progression, r, horloge(), __VERSION__), null, 1))
      setMessage(`Sauvegarde « ${nom} » prête.`)
    } catch (e) {
      setMessage(`Export impossible : ${(e as Error).message}`)
    }
  }

  const importer = async (f: File) => {
    try {
      const lu = lireSauvegarde(await f.text())
      const n = Object.keys(lu.progression.volets).length
      if (!confirm(`Remplacer la progression de ce téléphone par celle du fichier (${n} volets révisés) ? L’actuelle sera perdue.`)) return
      remplacer(lu.progression, lu.reglages)
      setMessage('Progression importée.')
    } catch (e) {
      setMessage((e as Error).message)
    }
  }

  return (
    <main className="ecran">
      <h1 className="titre-ecran">Réglages</h1>
      {message && (
        <div className="bandeau info" role="status">
          {message}
        </div>
      )}

      <span className="etiquette">Niveaux</span>
      <div className="bloc">
        <label className="champ">
          <span>Niveau de départ</span>
          <select value={r.niveauDepart} onChange={(e) => maj({ niveauDepart: e.target.value as Niveau })}>
            {NIVEAUX.map((n) => (
              <option key={n} value={n}>
                {NOMS_NIVEAUX[n]}
              </option>
            ))}
          </select>
        </label>
        <label className="champ">
          <span>Niveaux en dessous du niveau de départ</span>
          <select value={r.niveauxInferieurs} onChange={(e) => maj({ niveauxInferieurs: e.target.value as TypeReglages['niveauxInferieurs'] })}>
            <option value="acquis">considérés comme acquis</option>
            <option value="masques">masqués</option>
          </select>
        </label>
        <span>Niveaux actifs (cartes neuves proposées)</span>
        <div className="puces" style={{ marginTop: 6 }}>
          {NIVEAUX.map((n) => (
            <button
              key={n}
              className="puce"
              aria-pressed={r.niveauxActifs.includes(n)}
              onClick={() => {
                const l = basculerListe(r.niveauxActifs, n)
                if (l.length) maj({ niveauxActifs: NIVEAUX.filter((x) => l.includes(x)) })
              }}
            >
              {NOMS_NIVEAUX[n]}
            </button>
          ))}
        </div>
      </div>

      <span className="etiquette">Affichage</span>
      <div className="bloc">
        <label className="champ">
          <span>Thème</span>
          <select value={r.theme} onChange={(e) => maj({ theme: e.target.value as TypeReglages['theme'] })}>
            <option value="auto">comme le téléphone</option>
            <option value="clair">clair</option>
            <option value="sombre">sombre</option>
          </select>
        </label>
      </div>

      <span className="etiquette">Sauvegarde</span>
      <div className="bloc">
        <p className="petit">
          La progression n’existe que sur ce téléphone : l’export en JSON est la seule sauvegarde. Il contient aussi les
          réglages.
        </p>
        <button className="bouton" onClick={() => void exporter()}>
          Exporter la progression
        </button>
        <button className="bouton" onClick={() => fichier.current?.click()}>
          Importer une sauvegarde
        </button>
        <input
          ref={fichier}
          type="file"
          accept="application/json,.json"
          hidden
          onChange={(e) => {
            const f = e.target.files?.[0]
            if (f) void importer(f)
            e.target.value = ''
          }}
        />
      </div>

      <span className="etiquette">À propos</span>
      <div className="bloc petit">
        <p>
          {contenu.cartes.length} cartes, {contenu.exercices.length} exercices. Version {__VERSION__} ({__BUILD__}).
        </p>
        <p className="doux">
          Le découpage Terminale → L1 → L2 → L3 → agrégation est un choix de progression, pas un programme officiel (seuls la
          Terminale et l’agrégation en ont un).
        </p>
      </div>
    </main>
  )
}
