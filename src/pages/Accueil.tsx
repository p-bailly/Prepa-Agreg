// Accueil : ce qu'il y a à faire aujourd'hui, le gros bouton « Réviser », la série de jours et le niveau atteint.
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { resumeDuJour } from '../repetition/seance'
import { joursAvantEcrit, niveauAtteint, progresNiveaux, serieJours } from '../stats/tableau'

export function Accueil({ surReviser }: { surReviser: () => void }) {
  const { contenu, progression, reglages, horloge, aujourdhui } = useDonnees()
  const jour = aujourdhui()
  const r = resumeDuJour(contenu, progression, reglages, horloge(), jour)
  const serie = serieJours(progression, jour)
  const niveau = niveauAtteint(progresNiveaux(contenu, progression, reglages))
  const faites = progression.jours[jour]?.cartes ?? 0
  const rebours = joursAvantEcrit(reglages.dateEcrit, jour)
  const f = reglages.filtres
  const filtre = f.themes.length + f.niveaux.length + f.etats.length > 0
  const rien = r.aRevoir === 0 && r.nouvelles === 0
  return (
    <>
      <main className="ecran">
        <h1 className="titre-ecran">Prépa Agreg</h1>
        <span className="etiquette">Aujourd’hui</span>
        <div className="bloc">
          <p style={{ fontSize: 22, margin: 0 }}>
            <strong className="grand-nombre" style={{ fontSize: 40 }}>
              {r.aRevoir}
            </strong>{' '}
            {r.aRevoir > 1 ? 'cartes à revoir' : 'carte à revoir'}
          </p>
          <p style={{ fontSize: 20 }}>
            et {r.nouvelles} {r.nouvelles > 1 ? 'nouvelles' : 'nouvelle'}
          </p>
          <p className="petit doux">
            {faites} / {reglages.objectifQuotidien} cartes revues aujourd’hui (objectif quotidien)
            {filtre && ' · séance filtrée (voir Parcours)'}
          </p>
        </div>
        {rebours !== null && (
          <div className="bloc">
            <span className="etiquette">Écrit</span>
            <p style={{ margin: 0 }}>
              {rebours > 0 ? (
                <>
                  <strong className="grand-nombre" style={{ fontSize: 32 }}>
                    {rebours}
                  </strong>{' '}
                  {rebours > 1 ? 'jours' : 'jour'} avant l’écrit
                </>
              ) : rebours === 0 ? (
                'C’est aujourd’hui. Bon courage.'
              ) : (
                'La date de l’écrit est passée (Réglages pour la changer).'
              )}
            </p>
          </div>
        )}
        <div className="deux-colonnes">
          <div className="bloc">
            <span className="etiquette">Série</span>
            <div className="grand-nombre">{serie}</div>
            <div className="petit">{serie > 1 ? 'jours de suite' : 'jour'}</div>
          </div>
          <div className="bloc">
            <span className="etiquette">Niveau atteint</span>
            <div className="grand-nombre" style={{ fontSize: niveau && niveau.length > 3 ? 26 : 44 }}>
              {niveau ? NOMS_NIVEAUX[niveau] : '—'}
            </div>
            <div className="petit">80 % des cartes acquises</div>
          </div>
        </div>
      </main>
      <div className="pied">
        <div className="pied-interieur">
          <button className="bouton-principal" onClick={surReviser} disabled={rien}>
            Réviser
            <span style={{ fontSize: 15, letterSpacing: '0.02em', textTransform: 'none', fontWeight: 400 }}>
              {rien ? 'rien à revoir pour l’instant' : `environ ${r.minutes} minute${r.minutes > 1 ? 's' : ''}`}
            </span>
          </button>
        </div>
      </div>
    </>
  )
}
