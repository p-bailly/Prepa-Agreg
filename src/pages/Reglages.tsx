// Réglages : niveau de départ, niveaux actifs, modes de révision, séance et plafond de cartes neuves, objectif, date de
// l'écrit (facultative), thème, rappel quotidien, export et import de la progression, signalements.
import { useRef, useState, type ReactNode } from 'react'
import { NIVEAUX, NOMS_NIVEAUX, type Niveau } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { partagerFichier } from '../lib/fichiers'
import { annulerRappel, natif, programmerRappel } from '../lib/notifications'
import { MODES, NOMS_MODES, type Mode, type Reglages as TypeReglages } from '../lib/reglages'
import { fabriquerSauvegarde, lireSauvegarde } from '../lib/sauvegarde'

function Interrupteur({ actif, surChange, children }: { actif: boolean; surChange: (v: boolean) => void; children: ReactNode }) {
  return (
    <button className="interrupteur" role="switch" aria-checked={actif} onClick={() => surChange(!actif)}>
      <span>{children}</span>
      <span className="etat">{actif ? 'oui' : 'non'}</span>
    </button>
  )
}

function Nombre({ libelle, valeur, min, max, surChange }: { libelle: string; valeur: number; min: number; max: number; surChange: (n: number) => void }) {
  return (
    <label className="champ">
      <span>{libelle}</span>
      <input
        type="number"
        inputMode="numeric"
        min={min}
        max={max}
        value={valeur}
        onChange={(e) => {
          const n = Number(e.target.value)
          if (Number.isFinite(n)) surChange(Math.min(max, Math.max(min, Math.round(n))))
        }}
      />
    </label>
  )
}

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

  const changerRappel = async (actif: boolean, heure = r.rappel.heure) => {
    if (actif) {
      const ok = await programmerRappel(heure)
      if (!ok) {
        setMessage(natif ? 'Notifications refusées : autorisez-les dans les paramètres Android.' : 'Le rappel ne marche que dans l’appli Android.')
        maj({ rappel: { actif: false, heure } })
        return
      }
    } else await annulerRappel()
    maj({ rappel: { actif, heure } })
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

      <span className="etiquette">Modes de révision</span>
      <div className="bloc">
        {MODES.map((m) => (
          <Interrupteur
            key={m}
            actif={r.modesActifs.includes(m)}
            surChange={() => {
              const l = basculerListe<Mode>(r.modesActifs, m)
              if (l.length) maj({ modesActifs: MODES.filter((x) => l.includes(x)) })
            }}
          >
            {NOMS_MODES[m]}
          </Interrupteur>
        ))}
        <p className="petit doux" style={{ marginTop: 8 }}>
          Un mode n’est proposé que pour les cartes qui ont ce qu’il faut ; les modes alternent d’une carte à l’autre.
        </p>
      </div>

      <span className="etiquette">Séance</span>
      <div className="bloc">
        <Nombre libelle="Cartes par séance (au plus)" valeur={r.seanceCartes} min={1} max={100} surChange={(n) => maj({ seanceCartes: n })} />
        <Nombre libelle="Durée d’une séance (minutes, au plus)" valeur={r.seanceMinutes} min={1} max={120} surChange={(n) => maj({ seanceMinutes: n })} />
        <Nombre libelle="Cartes neuves par jour (plafond)" valeur={r.plafondNeuves} min={0} max={100} surChange={(n) => maj({ plafondNeuves: n })} />
        <Nombre libelle="Objectif quotidien (cartes revues)" valeur={r.objectifQuotidien} min={1} max={500} surChange={(n) => maj({ objectifQuotidien: n })} />
      </div>

      <span className="etiquette">Écrit</span>
      <div className="bloc">
        <label className="champ">
          <span>Date de l’écrit (facultative ; sans date, pas de compte à rebours)</span>
          <input type="date" value={r.dateEcrit ?? ''} onChange={(e) => maj({ dateEcrit: e.target.value || null })} />
        </label>
        {r.dateEcrit && (
          <button className="bouton compact" onClick={() => maj({ dateEcrit: null })}>
            Effacer la date
          </button>
        )}
      </div>

      <span className="etiquette">Affichage et rappel</span>
      <div className="bloc">
        <label className="champ">
          <span>Thème</span>
          <select value={r.theme} onChange={(e) => maj({ theme: e.target.value as TypeReglages['theme'] })}>
            <option value="auto">comme le téléphone</option>
            <option value="clair">clair</option>
            <option value="sombre">sombre</option>
          </select>
        </label>
        <Interrupteur actif={r.rappel.actif} surChange={(v) => void changerRappel(v)}>
          Rappel quotidien
        </Interrupteur>
        {r.rappel.actif && (
          <label className="champ" style={{ marginTop: 10 }}>
            <span>Heure du rappel</span>
            <input type="time" value={r.rappel.heure} onChange={(e) => e.target.value && void changerRappel(true, e.target.value)} />
          </label>
        )}
      </div>

      <span className="etiquette">Sauvegarde</span>
      <div className="bloc">
        <p className="petit">
          La progression n’existe que sur ce téléphone : l’export en JSON est la seule sauvegarde. Il contient aussi les
          réglages et les erreurs signalées ({progression.signalements.length}).
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

      {progression.signalements.length > 0 && (
        <>
          <span className="etiquette">Erreurs signalées</span>
          <div className="bloc petit">
            <ul style={{ margin: 0, paddingLeft: '1.1em' }}>
              {progression.signalements.map((s, i) => (
                <li key={i}>
                  <code>{s.cible}</code>
                  {s.ou ? ` (${s.ou})` : ''} : {s.commentaire}
                </li>
              ))}
            </ul>
          </div>
        </>
      )}

      <span className="etiquette">À propos</span>
      <div className="bloc petit">
        <p>
          <strong>Contenu rédigé et contrôlé par une IA, non relu par un humain.</strong> En cas de doute, une carte affiche
          « Non garanti » ; signalez toute erreur depuis la carte ou l’exercice.
        </p>
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
