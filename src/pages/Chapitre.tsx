// Cartes d'un chapitre, avec leur état ; toucher une carte l'ouvre en révision libre.
import { NOMS_NIVEAUX } from '../contenu/types'
import { Pastille } from '../composants/Icones'
import { useDonnees } from '../etat/Donnees'
import { etatAffiche, prerequisManquants, sousLeDepart } from '../repetition/etat'

const SYMBOLES = { neuve: 'vide', 'en-cours': 'moitie', acquise: 'plein' } as const
const NOMS = { neuve: 'neuve', 'en-cours': 'en cours', acquise: 'acquise' } as const

export function Chapitre({ id, surCarte, surRetour }: { id: string; surCarte: (id: string) => void; surRetour: () => void }) {
  const { contenu, progression, reglages, cartesParId } = useDonnees()
  const theme = contenu.themes.find((t) => t.chapitres.some((c) => c.id === id))
  const ch = theme?.chapitres.find((c) => c.id === id)
  const cartes = contenu.cartes.filter((k) => k.chapitre === id)
  return (
    <main className="ecran">
      <button className="retour" onClick={surRetour}>
        ← Parcours
      </button>
      <h1 className="titre-ecran">{ch?.titre}</h1>
      <span className="etiquette">
        {theme?.titre} · {ch ? NOMS_NIVEAUX[ch.niveau] : ''}
      </span>
      <ul className="liste-cartes">
        {cartes.map((k) => {
          const e = etatAffiche(k, progression.volets, reglages)
          const manquants = e === 'neuve' && !sousLeDepart(k, reglages) ? prerequisManquants(k, cartesParId, progression.volets, reglages) : []
          return (
            <li key={k.id}>
              <button className="bouton" style={{ justifyContent: 'flex-start', textAlign: 'left' }} onClick={() => surCarte(k.id)}>
                <Pastille etat={SYMBOLES[e]} />
                <span>
                  {k.titre}
                  <span className="petit doux" style={{ display: 'block' }}>
                    {NOMS[e]}
                    {manquants.length > 0 && ` · d’abord : ${manquants.map((m) => cartesParId.get(m)?.titre ?? m).join(', ')}`}
                  </span>
                </span>
              </button>
            </li>
          )
        })}
      </ul>
      <p className="petit doux">Une carte ouverte ici est révisée hors séance : la notation compte comme une révision.</p>
    </main>
  )
}
