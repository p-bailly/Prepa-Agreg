// Lignes des listes (chapitres, notions, exercices) : pastille d'avancement, intitulé (formules permises), ligne grise
// dessous. Toute la ligne est un bouton.
import type { ReactNode } from 'react'
import type { Avancement } from '../chapitres/plan'
import type { Exercice } from '../contenu/types'
import type { EtatExercice } from '../exercices/niveaux'
import { Pastille } from './Icones'
import { Texte } from './Texte'

const PASTILLES = { vide: 'vide', 'a-faire': 'vide', 'en-cours': 'moitie', fait: 'plein' } as const
export const NOMS_AVANCEMENT: Record<Avancement, string> = { vide: 'pas encore rédigé', 'a-faire': 'à faire', 'en-cours': 'en cours', fait: 'fait' }

export function Ligne({
  avancement,
  numero,
  titre,
  dessous,
  inactif,
  surClic,
}: {
  avancement: Avancement
  numero?: string
  titre: string
  dessous: ReactNode
  inactif?: boolean
  surClic: () => void
}) {
  return (
    <button className="ligne-liste" disabled={inactif} onClick={surClic}>
      <span className="pastille" aria-hidden="true">
        <Pastille etat={PASTILLES[avancement]} />
      </span>
      {numero && <span className="numero">{numero}</span>}
      <span className="intitule">
        <Texte texte={titre} enLigne />
        <small>{dessous}</small>
      </span>
    </button>
  )
}

const AVANCEMENT_EXERCICE: Record<EtatExercice, Avancement> = { 'a-faire': 'a-faire', partiel: 'en-cours', reussi: 'fait' }
const NOMS_EXERCICE: Record<EtatExercice, string> = { 'a-faire': 'à faire', partiel: 'tenté', reussi: 'réussi' }

export function LigneExercice({ exercice, etat, avant, surExercice }: { exercice: Exercice; etat: EtatExercice; avant?: string; surExercice: (id: string) => void }) {
  return (
    <Ligne
      avancement={AVANCEMENT_EXERCICE[etat]}
      titre={exercice.titre}
      dessous={`${avant ? `${avant} · ` : ''}difficulté ${exercice.difficulte}/5 · ${NOMS_EXERCICE[etat]}`}
      surClic={() => surExercice(exercice.id)}
    />
  )
}
