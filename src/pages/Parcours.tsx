// Parcours : choix d'un thème, puis chemin vertical de ses chapitres, du plus bas niveau au plus haut, avec leur numéro et
// leur avancement (● fait, ◐ en cours, ○ à faire ou pas encore rédigé). Toucher un chapitre rédigé l'ouvre.
import { useMemo, useState } from 'react'
import { avancementChapitre, chapitresNumerotes } from '../chapitres/plan'
import { NOMS_AVANCEMENT } from '../composants/Lignes'
import { Pastille } from '../composants/Icones'
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'

const PASTILLES = { vide: 'vide', 'a-faire': 'vide', 'en-cours': 'moitie', fait: 'plein' } as const

export function Parcours({ surChapitre }: { surChapitre: (id: string) => void }) {
  const { contenu, progression } = useDonnees()
  const [themeId, setThemeId] = useState(() => contenu.themes.find((t) => contenu.cartes.some((k) => k.theme === t.id))?.id ?? contenu.themes[0].id)
  const chapitres = useMemo(() => chapitresNumerotes(contenu).filter((x) => x.theme.id === themeId), [contenu, themeId])

  return (
    <main className="ecran">
      <h1 className="titre-ecran">Parcours</h1>
      <label className="champ">
        <span className="etiquette">Thème</span>
        <select value={themeId} onChange={(e) => setThemeId(e.target.value)}>
          {contenu.themes.map((t) => {
            const n = contenu.cartes.filter((k) => k.theme === t.id).length
            return (
              <option key={t.id} value={t.id}>
                {t.titre} {n ? `(${n} notions)` : '(à venir)'}
              </option>
            )
          })}
        </select>
      </label>

      <ol className="chemin">
        {chapitres.map(({ numero, chapitre }, i) => {
          const a = avancementChapitre(contenu, chapitre.id, progression)
          const nouveauNiveau = i === 0 || chapitres[i - 1].chapitre.niveau !== chapitre.niveau
          return (
            <li key={chapitre.id}>
              {nouveauNiveau && <div className="separateur-niveau">{NOMS_NIVEAUX[chapitre.niveau]}</div>}
              <div className="noeud">
                <span className="pastille" aria-hidden="true">
                  <Pastille etat={PASTILLES[a]} />
                </span>
                <button
                  className="intitule"
                  style={{ background: 'none', border: 'none', textAlign: 'left', padding: '6px 0', minHeight: 44 }}
                  disabled={a === 'vide'}
                  onClick={() => surChapitre(chapitre.id)}
                >
                  {numero}. {chapitre.titre}
                  <small>{a === 'vide' ? `${NOMS_AVANCEMENT.vide} (${chapitre.prevues} notions prévues)` : NOMS_AVANCEMENT[a]}</small>
                </button>
              </div>
            </li>
          )
        })}
      </ol>
    </main>
  )
}
