// Les quatre boutons de notation d'un volet, avec le délai avant la prochaine révision sous chacun ; chaque note a sa
// couleur (classe `note-n`), toujours doublée de son nom.
import { delaiAffiche, NOMS_NOTES, NOTES, type EtatVolet, type Note } from '../repetition/planification'

interface Props {
  etat: EtatVolet | undefined
  maintenant: number
  jour: string
  surNote: (n: Note) => void
}

export function BoutonsNotation({ etat, maintenant, jour, surNote }: Props) {
  return (
    <div className="rangee-boutons notation">
      {NOTES.map((n) => (
        <button key={n} className={`bouton note-${n}`} onClick={() => surNote(n)}>
          <span>{NOMS_NOTES[n]}</span>
          <span className="delai">{delaiAffiche(etat, n, maintenant, jour)}</span>
        </button>
      ))}
    </div>
  )
}
