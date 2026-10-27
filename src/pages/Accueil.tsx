// Accueil : ce qu'il y a à faire aujourd'hui avec le bouton « Réviser » dans le même bloc, puis la série de jours, le
// niveau atteint, le compte à rebours de l'écrit et l'objectif quotidien.
import { Lecture } from '../composants/Icones'
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { resumeDuJour } from '../repetition/seance'
import { joursAvantEcrit, niveauAtteint, progresNiveaux, serieJours } from '../stats/tableau'

/** « mar. 27 oct. » pour un jour civil `AAAA-MM-JJ` (midi UTC : le jour ne bascule dans aucun fuseau). */
const dateCourte = (jour: string) =>
  new Date(`${jour}T12:00:00Z`).toLocaleDateString('fr-FR', { weekday: 'short', day: 'numeric', month: 'short', timeZone: 'UTC' })

export function Accueil({ surReviser }: { surReviser: () => void }) {
  const { contenu, progression, reglages, horloge, aujourdhui } = useDonnees()
  const jour = aujourdhui()
  const r = resumeDuJour(contenu, progression, reglages, horloge(), jour)
  const serie = serieJours(progression, jour)
  const niveau = niveauAtteint(progresNiveaux(contenu, progression, reglages))
  const nomNiveau = niveau ? NOMS_NIVEAUX[niveau] : '—'
  const faites = progression.jours[jour]?.cartes ?? 0
  const rebours = joursAvantEcrit(reglages.dateEcrit, jour)
  const f = reglages.filtres
  const filtre = f.themes.length + f.niveaux.length + f.etats.length > 0
  const rien = r.aRevoir === 0 && r.nouvelles === 0
  const objectif = Math.max(1, reglages.objectifQuotidien)
  return (
    <main className="ecran">
      <div className="entete-accueil">
        <h1 className="titre-ecran">Prépa Agreg</h1>
        <span className="petit doux">{dateCourte(jour)}</span>
      </div>

      <section className="hero" aria-label="Aujourd’hui">
        <span className="hero-etiquette">Aujourd’hui</span>
        <div className="hero-compte">
          <strong>{r.aRevoir}</strong>
          <span>
            {r.aRevoir > 1 ? 'cartes à revoir' : 'carte à revoir'}
            <br />
            <span className="hero-doux">
              + {r.nouvelles} {r.nouvelles > 1 ? 'nouvelles' : 'nouvelle'}
            </span>
          </span>
        </div>
        <button className="bouton-hero" onClick={surReviser} disabled={rien}>
          Réviser
          <span className="duree">
            {rien ? 'rien à revoir pour l’instant' : `environ ${r.minutes} min`}
            <span className="rond">
              <Lecture />
            </span>
          </span>
        </button>
      </section>

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

      <div className="bloc">
        <div className="ligne" style={{ marginBottom: 8 }}>
          <span className="petit" style={{ fontWeight: 700 }}>
            Objectif du jour
          </span>
          <span className="petit" style={{ fontWeight: 700 }}>
            {faites} / {reglages.objectifQuotidien}
          </span>
        </div>
        <div className="barre" role="img" aria-label={`${faites} cartes revues sur ${reglages.objectifQuotidien}`}>
          <span style={{ width: `${Math.min(100, (100 * faites) / objectif)}%` }} />
        </div>
        {filtre && <p className="petit doux" style={{ marginTop: 10 }}>Séance filtrée (voir Parcours).</p>}
        {rebours !== null && rebours < 0 && (
          <p className="petit doux" style={{ marginTop: 10 }}>
            La date de l’écrit est passée (Réglages pour la changer).
          </p>
        )}
      </div>
    </main>
  )
}
