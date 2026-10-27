// Accueil : ce qu'il y a à faire aujourd'hui (cartes à revoir, leçon du jour) avec le bouton principal dans le même bloc
// (« Réviser », ou « Apprendre » quand rien n'est à revoir : on apprend avant de réviser), puis la leçon du jour, la série
// de jours, le niveau atteint, le compte à rebours de l'écrit et l'objectif quotidien.
import { Lecture } from '../composants/Icones'
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { resumeDuJour, type Lecon } from '../repetition/seance'
import { joursAvantEcrit, niveauAtteint, progresNiveaux, serieJours } from '../stats/tableau'

/** « mar. 27 oct. » pour un jour civil `AAAA-MM-JJ` (midi UTC : le jour ne bascule dans aucun fuseau). */
const dateCourte = (jour: string) =>
  new Date(`${jour}T12:00:00Z`).toLocaleDateString('fr-FR', { weekday: 'short', day: 'numeric', month: 'short', timeZone: 'UTC' })

export function Accueil({ surReviser, surApprendre }: { surReviser: () => void; surApprendre: (l: Lecon) => void }) {
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
  const objectif = Math.max(1, reglages.objectifQuotidien)
  const lecon = r.lecon
  const theme = lecon && contenu.themes.find((t) => t.chapitres.some((c) => c.id === lecon.chapitre.id))
  // On apprend d'abord : sans révision due, le bouton principal lance la leçon du jour.
  const principal: 'reviser' | 'apprendre' | null = r.aRevoir > 0 ? 'reviser' : lecon ? 'apprendre' : null
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
              + {r.aApprendre} {r.aApprendre > 1 ? 'notions à apprendre' : 'notion à apprendre'}
            </span>
          </span>
        </div>
        {principal === 'apprendre' && lecon ? (
          <button className="bouton-hero" onClick={() => surApprendre(lecon)}>
            Apprendre
            <span className="duree">
              {`leçon du jour, environ ${r.minutesLecon} min`}
              <span className="rond">
                <Lecture />
              </span>
            </span>
          </button>
        ) : (
          <button className="bouton-hero" onClick={surReviser} disabled={principal === null}>
            Réviser
            <span className="duree">
              {principal === null ? 'rien à revoir ni à apprendre' : `environ ${r.minutesRevision} min`}
              <span className="rond">
                <Lecture />
              </span>
            </span>
          </button>
        )}
      </section>

      {lecon && (
        <div className="bloc">
          <span className="etiquette">Leçon du jour</span>
          <p style={{ margin: '0 0 4px', fontWeight: 700 }}>{lecon.chapitre.titre}</p>
          <p className="petit doux">
            {theme?.titre} · {NOMS_NIVEAUX[lecon.chapitre.niveau]} · {lecon.cartes.length} notion
            {lecon.cartes.length > 1 ? 's' : ''} · environ {r.minutesLecon} min
          </p>
          {principal !== 'apprendre' && (
            <button className="bouton" onClick={() => surApprendre(lecon)}>
              Apprendre
            </button>
          )}
        </div>
      )}

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
        {filtre && <p className="petit doux" style={{ marginTop: 10 }}>Filtres actifs (voir Parcours).</p>}
        {rebours !== null && rebours < 0 && (
          <p className="petit doux" style={{ marginTop: 10 }}>
            La date de l’écrit est passée (Réglages pour la changer).
          </p>
        )}
      </div>
    </main>
  )
}
