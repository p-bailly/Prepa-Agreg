// Parcours : choix d'un thème, cinq jauges de niveau, filtres de séance, puis chemin vertical des chapitres du plus bas
// niveau au plus haut (● acquis, ◐ en cours, ○ à venir, cadenas si les prérequis manquent). Toucher un nœud ouvre ses cartes.
import { useState } from 'react'
import { Volet } from '../composants/Blocs'
import { Cadenas, Pastille } from '../composants/Icones'
import { NIVEAUX, NOMS_NIVEAUX, rangNiveau, type Carte, type Chapitre, type Contenu } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import type { Progression } from '../lib/progression'
import { ETATS_CARTE, NOMS_ETATS, type Reglages } from '../lib/reglages'
import { etatAffiche, prerequisManquants, sousLeDepart } from '../repetition/etat'

export type EtatChapitre = 'acquis' | 'en-cours' | 'a-venir' | 'bloque' | 'vide' | 'masque'

export function etatChapitre(
  ch: Chapitre,
  c: Contenu,
  p: Progression,
  r: Reglages,
  parId: Map<string, Carte>,
): { etat: EtatChapitre; manquants: string[]; cartes: Carte[] } {
  const cartes = c.cartes.filter((k) => k.chapitre === ch.id)
  if (rangNiveau(ch.niveau) < rangNiveau(r.niveauDepart) && r.niveauxInferieurs === 'masques')
    return { etat: 'masque', manquants: [], cartes }
  if (!cartes.length) return { etat: 'vide', manquants: [], cartes }
  const etats = cartes.map((k) => etatAffiche(k, p, r))
  if (etats.every((e) => e === 'acquise')) return { etat: 'acquis', manquants: [], cartes }
  if (etats.some((e) => e !== 'neuve')) return { etat: 'en-cours', manquants: [], cartes }
  const manquants = cartes.map((k) => (sousLeDepart(k, r) ? [] : prerequisManquants(k, parId, p, r)))
  if (manquants.every((m) => m.length > 0)) return { etat: 'bloque', manquants: [...new Set(manquants.flat())], cartes }
  return { etat: 'a-venir', manquants: [], cartes }
}

const PASTILLES = {
  acquis: 'plein',
  'en-cours': 'moitie',
  'a-venir': 'vide',
  vide: 'vide',
  masque: 'tiret',
} as const

export function Parcours({ surChapitre }: { surChapitre: (id: string) => void }) {
  const { contenu, progression, reglages, changerReglages, cartesParId } = useDonnees()
  const [themeId, setThemeId] = useState(() => contenu.themes.find((t) => contenu.cartes.some((k) => k.theme === t.id))?.id ?? contenu.themes[0].id)
  const [cadenas, setCadenas] = useState<string | null>(null)
  const theme = contenu.themes.find((t) => t.id === themeId)!
  const chapitres = [...theme.chapitres].sort((a, b) => rangNiveau(a.niveau) - rangNiveau(b.niveau))
  const f = reglages.filtres

  const basculer = <K extends keyof Reglages['filtres']>(cle: K, v: Reglages['filtres'][K][number]) =>
    changerReglages((r) => {
      const l = r.filtres[cle] as string[]
      const n = l.includes(v) ? l.filter((x) => x !== v) : [...l, v]
      return { ...r, filtres: { ...r.filtres, [cle]: n } }
    })
  const nbFiltres = f.themes.length + f.niveaux.length + f.etats.length

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
                {t.titre} {n ? `(${n} cartes)` : '(à venir)'}
              </option>
            )
          })}
        </select>
      </label>

      <div className="bloc">
        {NIVEAUX.map((n) => {
          const visees = theme.chapitres.filter((ch) => ch.niveau === n).reduce((s, ch) => s + ch.prevues, 0)
          const cartes = contenu.cartes.filter((k) => k.theme === themeId && k.niveau === n)
          const acquises = cartes.filter((k) => etatAffiche(k, progression, reglages) === 'acquise').length
          const total = Math.max(visees, cartes.length)
          return (
            <div key={n} className="jauge">
              <div className="ligne">
                <span>{NOMS_NIVEAUX[n]}</span>
                <span className="petit">{total ? `${acquises} / ${total}` : 'rien de prévu'}</span>
              </div>
              <div className="barre">
                <span style={{ width: total ? `${(100 * acquises) / total}%` : 0 }} />
              </div>
            </div>
          )
        })}
      </div>

      <Volet titre="Filtres de séance" marque={nbFiltres ? `${nbFiltres} actif${nbFiltres > 1 ? 's' : ''}` : undefined}>
        <p className="petit doux" style={{ marginTop: 10 }}>
          La séance du jour ne prend que les cartes qui passent ces filtres (aucun choix = pas de filtre).
        </p>
        <span className="etiquette">Thèmes</span>
        <div className="puces">
          {contenu.themes.map((t) => (
            <button key={t.id} className="puce" aria-pressed={f.themes.includes(t.id)} onClick={() => basculer('themes', t.id)}>
              {t.titre}
            </button>
          ))}
        </div>
        <span className="etiquette">Niveaux</span>
        <div className="puces">
          {NIVEAUX.map((n) => (
            <button key={n} className="puce" aria-pressed={f.niveaux.includes(n)} onClick={() => basculer('niveaux', n)}>
              {NOMS_NIVEAUX[n]}
            </button>
          ))}
        </div>
        <span className="etiquette">États</span>
        <div className="puces">
          {ETATS_CARTE.map((e) => (
            <button key={e} className="puce" aria-pressed={f.etats.includes(e)} onClick={() => basculer('etats', e)}>
              {NOMS_ETATS[e]}
            </button>
          ))}
        </div>
        {nbFiltres > 0 && (
          <button className="bouton compact" onClick={() => changerReglages((r) => ({ ...r, filtres: { themes: [], niveaux: [], etats: [] } }))}>
            Retirer les filtres
          </button>
        )}
      </Volet>

      <ol className="chemin">
        {chapitres.map((ch, i) => {
          const e = etatChapitre(ch, contenu, progression, reglages, cartesParId)
          const nouveauNiveau = i === 0 || chapitres[i - 1].niveau !== ch.niveau
          const libelle =
            e.etat === 'vide'
              ? `à venir — pas encore rédigé (${ch.prevues} cartes prévues)`
              : e.etat === 'masque'
                ? 'masqué (sous le niveau de départ)'
                : `${e.cartes.length} carte${e.cartes.length > 1 ? 's' : ''}${e.etat === 'acquis' ? ' · acquis' : e.etat === 'en-cours' ? ' · en cours' : e.etat === 'bloque' ? ' · prérequis manquants' : ''}`
          return (
            <li key={ch.id}>
              {nouveauNiveau && <div className="separateur-niveau">{NOMS_NIVEAUX[ch.niveau]}</div>}
              <div className="noeud">
                {e.etat === 'bloque' ? (
                  <button
                    className="pastille"
                    aria-label="Prérequis manquants : voir lesquels"
                    onClick={() => setCadenas(cadenas === ch.id ? null : ch.id)}
                  >
                    <Cadenas />
                  </button>
                ) : (
                  <span className="pastille" aria-hidden="true">
                    <Pastille etat={PASTILLES[e.etat]} />
                  </span>
                )}
                <button
                  className="intitule"
                  style={{ background: 'none', border: 'none', textAlign: 'left', padding: '6px 0', minHeight: 44 }}
                  disabled={e.etat === 'vide' || e.etat === 'masque'}
                  onClick={() => surChapitre(ch.id)}
                >
                  {ch.titre} ({NOMS_NIVEAUX[ch.niveau]})
                  <small>{libelle}</small>
                </button>
              </div>
              {cadenas === ch.id && (
                <div className="bloc petit" style={{ marginLeft: 58 }}>
                  À voir d’abord (au moins « en cours ») :{' '}
                  {e.manquants.map((id) => cartesParId.get(id)?.titre ?? id).join(' ; ')}.
                </div>
              )}
            </li>
          )
        })}
      </ol>
    </main>
  )
}
