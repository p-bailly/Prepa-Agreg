// Accueil : série de jours, niveau atteint, compte à rebours de l'écrit, puis la pyramide des chapitres (Terminale à la
// base, agrégation au sommet, reliés par leurs prérequis), qui mène à chaque chapitre.
import { Pyramide } from '../composants/Pyramide'
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { joursAvantEcrit, niveauAtteint, progresNiveaux, serieJours } from '../stats/tableau'

/** « mar. 27 oct. » pour un jour civil `AAAA-MM-JJ` (midi UTC : le jour ne bascule dans aucun fuseau). */
const dateCourte = (jour: string) =>
  new Date(`${jour}T12:00:00Z`).toLocaleDateString('fr-FR', { weekday: 'short', day: 'numeric', month: 'short', timeZone: 'UTC' })

export function Accueil({ surChapitre }: { surChapitre: (id: string) => void }) {
  const { contenu, progression, reglages, aujourdhui } = useDonnees()
  const jour = aujourdhui()
  const serie = serieJours(progression, jour)
  const niveau = niveauAtteint(progresNiveaux(contenu, progression, reglages))
  const nomNiveau = niveau ? NOMS_NIVEAUX[niveau] : '—'
  const rebours = joursAvantEcrit(reglages.dateEcrit, jour)
  return (
    <main className="ecran">
      <div className="entete-accueil">
        <h1 className="titre-ecran">Prépa Agreg</h1>
        <span className="petit doux">{dateCourte(jour)}</span>
      </div>

      <div className="tuiles">
        <div className="tuile sable">
          <span className="nom">Série</span>
          <span className="valeur">{serie} j</span>
        </div>
        <div className="tuile menthe">
          <span className="nom">Niveau atteint</span>
          <span className={nomNiveau.length > 3 ? 'valeur longue' : 'valeur'}>{nomNiveau}</span>
        </div>
        {rebours !== null && (
          <div className="tuile rose">
            <span className="nom">Écrit</span>
            <span className={rebours > 0 ? 'valeur' : 'valeur longue'}>
              {rebours > 0 ? `${rebours} j` : rebours === 0 ? 'aujourd’hui' : 'passé'}
            </span>
          </div>
        )}
      </div>

      <h2 className="sous-titre">Carte des notions</h2>
      <Pyramide surChapitre={surChapitre} />
    </main>
  )
}
