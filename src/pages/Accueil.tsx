// Accueil : ce qu'il y a à faire aujourd'hui (cartes à revoir, leçon du jour), le gros bouton (« Réviser », ou
// « Apprendre » quand rien n'est à revoir), la série de jours et le niveau atteint.
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { resumeDuJour, type Lecon } from '../repetition/seance'
import { joursAvantEcrit, niveauAtteint, progresNiveaux, serieJours } from '../stats/tableau'

export function Accueil({ surReviser, surApprendre }: { surReviser: () => void; surApprendre: (l: Lecon) => void }) {
  const { contenu, progression, reglages, horloge, aujourdhui } = useDonnees()
  const jour = aujourdhui()
  const r = resumeDuJour(contenu, progression, reglages, horloge(), jour)
  const serie = serieJours(progression, jour)
  const niveau = niveauAtteint(progresNiveaux(contenu, progression, reglages))
  const faites = progression.jours[jour]?.cartes ?? 0
  const rebours = joursAvantEcrit(reglages.dateEcrit, jour)
  const f = reglages.filtres
  const filtre = f.themes.length + f.niveaux.length + f.etats.length > 0
  const theme = r.lecon && contenu.themes.find((t) => t.chapitres.some((c) => c.id === r.lecon!.chapitre.id))
  // On apprend d'abord : sans révision due, le gros bouton lance la leçon.
  const principal: 'reviser' | 'apprendre' | null = r.aRevoir > 0 ? 'reviser' : r.lecon ? 'apprendre' : null
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
            et {r.aApprendre} {r.aApprendre > 1 ? 'notions nouvelles à apprendre' : 'notion nouvelle à apprendre'}
          </p>
          <p className="petit doux">
            {faites} / {reglages.objectifQuotidien} cartes revues aujourd’hui (objectif quotidien)
            {filtre && ' · filtres actifs (voir Parcours)'}
          </p>
        </div>
        {r.lecon && (
          <div className="bloc">
            <span className="etiquette">Leçon du jour</span>
            <p style={{ margin: '0 0 4px', fontSize: 20 }}>
              <strong>{r.lecon.chapitre.titre}</strong>
            </p>
            <p className="petit doux">
              {theme?.titre} · {NOMS_NIVEAUX[r.lecon.chapitre.niveau]} · {r.lecon.cartes.length} notion
              {r.lecon.cartes.length > 1 ? 's' : ''} · environ {r.minutesLecon} minutes
            </p>
            {principal !== 'apprendre' && (
              <button className="bouton" onClick={() => surApprendre(r.lecon!)}>
                Apprendre
              </button>
            )}
          </div>
        )}
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
          {principal === 'apprendre' ? (
            <button className="bouton-principal" onClick={() => surApprendre(r.lecon!)}>
              Apprendre
              <span style={{ fontSize: 15, letterSpacing: '0.02em', textTransform: 'none', fontWeight: 400 }}>
                leçon du jour, environ {r.minutesLecon} minutes
              </span>
            </button>
          ) : (
            <button className="bouton-principal" onClick={surReviser} disabled={principal === null}>
              Réviser
              <span style={{ fontSize: 15, letterSpacing: '0.02em', textTransform: 'none', fontWeight: 400 }}>
                {principal === null
                  ? 'rien à revoir ni à apprendre pour l’instant'
                  : `environ ${r.minutesRevision} minute${r.minutesRevision > 1 ? 's' : ''}`}
              </span>
            </button>
          )}
        </div>
      </div>
    </>
  )
}
