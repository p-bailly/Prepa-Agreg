// Leçon : on APPREND avant de réviser. Introduction du chapitre, puis une page par notion : texte de cours, énoncé
// encadré, exemple, démonstration, contre-exemples, pièges (vrai/faux commentés) et exercices d'application. Passer à
// la notion suivante la marque comme apprise : elle entre en révision espacée à partir du lendemain.
import { useState, type CSSProperties } from 'react'
import { BandeauDoute, Encadre, Volet } from '../composants/Blocs'
import { Signaler } from '../composants/Signalement'
import { Texte } from '../composants/Texte'
import { NOMS_NIVEAUX, type Carte, type Chapitre } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { apprendre } from '../lib/progression'
import { VoletExercices } from './EcranCarte'

interface Props {
  chapitre: Chapitre
  cartes: Carte[]
  surExercice: (id: string) => void
  surFin: () => void
}

export function EcranLecon({ chapitre, cartes, surExercice, surFin }: Props) {
  const { contenu, changerProgression, aujourdhui } = useDonnees()
  const intro = contenu.lecons?.find((l) => l.chapitre === chapitre.id)
  const theme = contenu.themes.find((t) => t.chapitres.some((c) => c.id === chapitre.id))
  // Pages : introduction (si elle existe), une par notion, puis le bilan.
  const debut = intro ? -1 : 0
  const [page, setPage] = useState(debut)
  const [apprises, setApprises] = useState<string[]>([])
  const carte = page >= 0 ? cartes[page] : undefined
  const fin = page >= cartes.length

  const aller = (n: number) => {
    // Quitter une page de notion vers la suite = notion apprise.
    if (carte && n > page && !apprises.includes(carte.id)) {
      changerProgression((p) => apprendre(p, [carte.id], aujourdhui()))
      setApprises((a) => [...a, carte.id])
    }
    setPage(n)
    window.scrollTo(0, 0)
  }

  const style = { '--pied': '100px' } as CSSProperties
  return (
    <div className="appli" style={style}>
      <main className="ecran">
        <div className="barre-haut">
          <button className="retour" onClick={surFin}>
            ← {fin ? 'Accueil' : 'Quitter'}
          </button>
          <span className="avancement">
            <span className="piste" aria-hidden="true">
              <span style={{ width: `${(100 * Math.min(cartes.length, Math.max(0, page + 1))) / Math.max(1, cartes.length)}%` }} />
            </span>
            {fin ? 'Bilan' : carte ? `${page + 1} / ${cartes.length}` : 'Introduction'}
          </span>
        </div>
        <header className="tete-carte">
          <span className="sur-titre">
            {fin ? 'Leçon terminée' : carte ? `Leçon · notion ${page + 1}` : 'Leçon'} · {theme?.titre} ·{' '}
            {NOMS_NIVEAUX[chapitre.niveau]}
          </span>
          <h1 className="titre-carte">
            <Texte texte={carte ? carte.titre : chapitre.titre} enLigne />
          </h1>
        </header>

        {!carte && !fin && (
          <>
            <BandeauDoute doute={intro?.doute} />
            {intro && <Texte texte={intro.introduction} />}
            <div className="bloc">
              <span className="etiquette">Au programme de cette leçon</span>
              <ol style={{ margin: 0, paddingLeft: '1.3em' }}>
                {cartes.map((c) => (
                  <li key={c.id}>
                    <Texte texte={c.titre} enLigne />
                  </li>
                ))}
              </ol>
            </div>
          </>
        )}

        {carte && (
          <>
            <BandeauDoute doute={carte.doute} />
            {carte.cours && <Texte texte={carte.cours} />}
            <Encadre genre={carte.type}>
              <Texte texte={carte.enonce} />
            </Encadre>
            {carte.exemple && (
              <Encadre genre="exemple">
                <Texte texte={carte.exemple} />
              </Encadre>
            )}
            {carte.demo && (
              <Volet titre="Démonstration">
                <div style={{ paddingTop: 12 }}>
                  <Texte texte={carte.demo} />
                </div>
              </Volet>
            )}
            {carte.contre_exemple && (
              <Volet titre="Pourquoi ces hypothèses">
                <Encadre genre="remarque" nom="Contre-exemples">
                  <Texte texte={carte.contre_exemple} />
                </Encadre>
              </Volet>
            )}
            {!!carte.affirmations?.length && (
              <Volet titre="Pièges à éviter">
                <div style={{ paddingTop: 12 }}>
                  {carte.affirmations.map((a, i) => (
                    <Encadre key={i} genre={a.vrai ? 'definition' : 'theoreme'} nom={a.vrai ? 'Vrai' : 'Faux'}>
                      <Texte texte={a.texte} />
                      {a.vrai && a.justification && <Texte texte={a.justification} />}
                      {!a.vrai && a.contre_exemple && (
                        <>
                          <p className="petit doux" style={{ marginTop: 8 }}>
                            Contre-exemple :
                          </p>
                          <Texte texte={a.contre_exemple} />
                        </>
                      )}
                    </Encadre>
                  ))}
                </div>
              </Volet>
            )}
            <Volet titre="S’exercer">
              <VoletExercices carte={carte} surExercice={surExercice} />
            </Volet>
            <Signaler cible={carte.id} genre="carte" ou="leçon" />
          </>
        )}

        {fin && (
          <>
            <div className="bloc">
              <p>
                <strong>{apprises.length}</strong> notion{apprises.length > 1 ? 's' : ''} apprise{apprises.length > 1 ? 's' : ''}.
              </p>
              <p className="doux">
                Elles entrent en révision espacée : première révision dès demain, dans la séance « Réviser ». La leçon reste
                accessible dans Parcours.
              </p>
            </div>
          </>
        )}
      </main>
      <div className="pied">
        <div className="pied-interieur">
          {fin ? (
            <button className="bouton-principal" onClick={surFin}>
              Terminer
            </button>
          ) : (
            <div style={{ display: 'grid', gridTemplateColumns: page > debut ? '1fr 2fr' : '1fr', gap: 8 }}>
              {page > debut && (
                <button className="bouton" style={{ minHeight: 76 }} onClick={() => aller(page - 1)}>
                  Précédente
                </button>
              )}
              <button className="bouton-principal" style={{ fontSize: 19 }} onClick={() => aller(page + 1)}>
                {!carte ? 'Commencer' : page === cartes.length - 1 ? 'Terminer' : 'Suivante'}
                {carte && (
                  <span style={{ fontSize: 14, letterSpacing: '0.02em', textTransform: 'none', fontWeight: 400 }}>
                    notion comprise
                  </span>
                )}
              </button>
            </div>
          )}
        </div>
      </div>
    </div>
  )
}
