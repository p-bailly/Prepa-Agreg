// Stats : série, cartes acquises, niveau atteint ; cartes revues par jour sur 30 jours ; réussite par thème (faibles en
// premier) ; exercices réussis par difficulté ; compte à rebours si une date d'écrit est saisie. Sans badge ni étoile.
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import {
  cartesAcquises,
  cartesParJour,
  exercicesParDifficulte,
  joursAvantEcrit,
  niveauAtteint,
  progresNiveaux,
  reussiteParTheme,
  serieJours,
} from '../stats/tableau'

const pourcent = (x: number) => `${Math.round(100 * x)} %`
/** Largeur CSS d'une barre (sans l'espace insécable du pourcentage affiché). */
const largeur = (x: number) => `${Math.round(1000 * Math.min(1, Math.max(0, x))) / 10}%`

export function Stats() {
  const { contenu, progression, reglages, aujourdhui } = useDonnees()
  const jour = aujourdhui()
  const niveaux = progresNiveaux(contenu, progression, reglages)
  const atteint = niveauAtteint(niveaux)
  const histo = cartesParJour(progression, jour, 30)
  const max = Math.max(1, ...histo.map((h) => h.cartes))
  const themes = reussiteParTheme(contenu, progression)
  const exercices = exercicesParDifficulte(contenu, progression)
  const rebours = joursAvantEcrit(reglages.dateEcrit, jour)
  const total = histo.reduce((s, h) => s + h.cartes, 0)
  return (
    <main className="ecran">
      <h1 className="titre-ecran">Stats</h1>
      <div className="trois-colonnes">
        <div>
          <div className="grand-nombre">{serieJours(progression, jour)}</div>
          <div className="petit">{serieJours(progression, jour) > 1 ? 'jours de suite' : 'jour de suite'}</div>
        </div>
        <div>
          <div className="grand-nombre">{cartesAcquises(contenu, progression, reglages)}</div>
          <div className="petit">cartes acquises</div>
        </div>
        <div>
          <div className="grand-nombre" style={{ fontSize: atteint && atteint.length > 3 ? 20 : 34 }}>
            {atteint ? NOMS_NIVEAUX[atteint] : '—'}
          </div>
          <div className="petit">niveau atteint</div>
        </div>
      </div>
      {rebours !== null && rebours >= 0 && (
        <p>
          <strong>{rebours}</strong> jour{rebours > 1 ? 's' : ''} avant l’écrit.
        </p>
      )}

      <span className="etiquette">Cartes revues par jour</span>
      <div className="bloc">
        <div className="histogramme" role="img" aria-label={`${total} cartes revues sur les 30 derniers jours`}>
          {histo.map((h, i) => (
            <span
              key={h.jour}
              className={i === histo.length - 1 ? 'aujourdhui' : undefined}
              style={{ height: `${(100 * h.cartes) / max}%` }}
              title={`${h.jour} : ${h.cartes}`}
            />
          ))}
        </div>
        <div className="ligne petit doux">
          <span>il y a 30 jours</span>
          <span>max {max === 1 && total === 0 ? 0 : max} / jour</span>
          <span>aujourd’hui</span>
        </div>
      </div>

      <span className="etiquette">Progression par niveau</span>
      <div className="bloc">
        {niveaux.map((n) => (
          <div key={n.niveau} className="jauge">
            <div className="ligne">
              <span>{NOMS_NIVEAUX[n.niveau]}</span>
              <span className="petit">
                {n.statut === 'acquis'
                  ? 'sous le niveau de départ : acquis'
                  : n.statut === 'masque'
                    ? 'masqué'
                    : `${n.acquises} / ${n.visees} (${n.existantes} rédigées)`}
              </span>
            </div>
            <div className="barre">
              <span style={{ width: largeur(n.part) }} />
            </div>
          </div>
        ))}
      </div>

      <span className="etiquette">Réussite par thème</span>
      <div className="bloc">
        {themes.length === 0 && <p className="doux">Pas encore de révision.</p>}
        {themes.map((t) => (
          <div key={t.theme} className="jauge">
            <div className="ligne">
              <span>{t.titre}</span>
              <span className="petit">
                {pourcent(t.part)} de {t.revisions}
              </span>
            </div>
            <div className="barre seconde">
              <span style={{ width: largeur(t.part) }} />
            </div>
          </div>
        ))}
        <p className="petit doux">Part des révisions notées Bien ou Facile.</p>
      </div>

      <span className="etiquette">Exercices réussis par difficulté</span>
      <div className="bloc">
        {exercices.map((e) => (
          <div key={e.difficulte} className="jauge">
            <div className="ligne">
              <span>Difficulté {e.difficulte}</span>
              <span className="petit">
                {e.reussis} / {e.total}
              </span>
            </div>
            <div className="barre">
              <span style={{ width: largeur(e.total ? e.reussis / e.total : 0) }} />
            </div>
          </div>
        ))}
      </div>
    </main>
  )
}
