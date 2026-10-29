// Cartes d'un chapitre, avec leur état ; toucher une carte l'ouvre en révision libre.
import { NOMS_NIVEAUX } from '../contenu/types'
import { Pastille } from '../composants/Icones'
import { Texte } from '../composants/Texte'
import { useDonnees } from '../etat/Donnees'
import { etatAffiche, prerequisManquants, sousLeDepart } from '../repetition/etat'

const SYMBOLES = { neuve: 'vide', 'en-cours': 'moitie', acquise: 'plein' } as const
const NOMS = { neuve: 'neuve', 'en-cours': 'en cours', acquise: 'acquise' } as const

export function Chapitre({
  id,
  surCarte,
  surLecon,
  surRetour,
}: {
  id: string
  surCarte: (id: string) => void
  surLecon: (chapitre: string, cartes: string[]) => void
  surRetour: () => void
}) {
  const { contenu, progression, reglages, cartesParId } = useDonnees()
  const theme = contenu.themes.find((t) => t.chapitres.some((c) => c.id === id))
  const ch = theme?.chapitres.find((c) => c.id === id)
  const cartes = contenu.cartes.filter((k) => k.chapitre === id)
  return (
    <main className="ecran">
      <button className="retour" onClick={surRetour}>
        ← Retour
      </button>
      <h1 className="titre-ecran">{ch?.titre}</h1>
      <span className="etiquette">
        {theme?.titre} · {ch ? NOMS_NIVEAUX[ch.niveau] : ''}
      </span>
      {cartes.length > 0 && (
        <button className="bouton" style={{ marginBottom: 14 }} onClick={() => surLecon(id, cartes.map((k) => k.id))}>
          {cartes.some((k) => etatAffiche(k, progression, reglages) === 'neuve') ? 'Apprendre la leçon' : 'Relire la leçon'}
        </button>
      )}
      <ul className="liste-cartes">
        {cartes.map((k) => {
          const e = etatAffiche(k, progression, reglages)
          const manquants = e === 'neuve' && !sousLeDepart(k, reglages) ? prerequisManquants(k, cartesParId, progression, reglages) : []
          return (
            <li key={k.id}>
              <button className="bouton" style={{ justifyContent: 'flex-start', textAlign: 'left' }} onClick={() => surCarte(k.id)}>
                <Pastille etat={SYMBOLES[e]} />
                <span>
                  <Texte texte={k.titre} enLigne />
                  <span className="petit doux" style={{ display: 'block' }}>
                    {NOMS[e]}
                    {manquants.length > 0 && (
                      <>
                        {' · d’abord : '}
                        <Texte texte={manquants.map((m) => cartesParId.get(m)?.titre ?? m).join(', ')} enLigne />
                      </>
                    )}
                  </span>
                </span>
              </button>
            </li>
          )
        })}
      </ul>
      <p className="petit doux">
        La leçon présente chaque notion avant qu’elle entre en révision. Une carte ouverte ici est révisée hors séance : la
        notation compte comme une révision.
      </p>
    </main>
  )
}
