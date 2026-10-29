// Icônes dessinées (pas d'emoji) : cadenas du parcours, onglets, thèmes, chevron des volets.
export function Cadenas({ taille = 20 }: { taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" strokeWidth="2">
      <rect x="5" y="11" width="14" height="10" />
      <path d="M8 11V7a4 4 0 0 1 8 0v4" />
    </svg>
  )
}

/** Pastille d'état : plein (acquis), moitié (en cours), vide (à venir), tiret (masqué). */
export function Pastille({ etat, taille = 22 }: { etat: 'plein' | 'moitie' | 'vide' | 'tiret'; taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" stroke="currentColor" strokeWidth="2">
      {etat === 'tiret' ? (
        <path d="M6 12h12" />
      ) : (
        <>
          <circle cx="12" cy="12" r="9" fill={etat === 'plein' ? 'currentColor' : 'none'} />
          {etat === 'moitie' && <path d="M12 3a9 9 0 0 1 0 18z" fill="currentColor" />}
        </>
      )}
    </svg>
  )
}

/** Case d'un niveau d'exercices : ■ réussi, ◧ fait en partie, □ à faire (dessinées : pas de glyphe à chercher). */
export function CaseNiveau({ etat, taille = 22 }: { etat: 'reussi' | 'partiel' | 'a-faire' | 'aucun'; taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" stroke="currentColor" strokeWidth="2">
      <rect x="3" y="3" width="18" height="18" fill={etat === 'reussi' ? 'currentColor' : 'none'} />
      {etat === 'partiel' && <rect x="3" y="3" width="9" height="18" fill="currentColor" />}
    </svg>
  )
}

const TRAITS = { fill: 'none', stroke: 'currentColor', strokeWidth: 1.9, strokeLinecap: 'round', strokeLinejoin: 'round' } as const

/** Icônes de la barre d'onglets. */
export function IconeOnglet({ onglet, taille = 22 }: { onglet: 'accueil' | 'parcours' | 'reglages'; taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" {...TRAITS}>
      {onglet === 'accueil' && <path d="M3 11l9-7 9 7v9h-6v-6H9v6H3z" />}
      {onglet === 'parcours' && <path d="M4 6h16M4 12h16M4 18h10" />}
      {onglet === 'reglages' && (
        <>
          <path d="M4 7h10M18 7h2M4 17h4M12 17h8" />
          <circle cx="16" cy="7" r="2" />
          <circle cx="10" cy="17" r="2" />
        </>
      )}
    </svg>
  )
}

/** Icône d'un thème du parcours (même identifiants que `content/themes.json`), tracée dans la couleur du thème. */
export function IconeTheme({ theme, taille = 18 }: { theme: string; taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" {...TRAITS} strokeWidth={2.1}>
      {theme === 'fondements' && (
        <>
          <circle cx="9" cy="12" r="5.5" />
          <circle cx="15" cy="12" r="5.5" />
        </>
      )}
      {theme === 'analyse-reelle' && <path d="M3 19c3 0 4-12 7-12s3 8 6 8 3-6 5-6" />}
      {theme === 'integration' && <path d="M15 4.5c-1.2-1.4-3.5-1-3.8 1.6l-1.4 11.8c-.3 2.6-2.6 3-3.8 1.6" />}
      {theme === 'equations-differentielles' && <path d="M3 12c2-5 4-5 6 0s4 5 6 0 4-5 6 0" />}
      {theme === 'algebre-lineaire' && (
        <>
          <path d="M7 4H5v16h2M17 4h2v16h-2" />
          <path d="M9.5 9h.01M14.5 9h.01M9.5 15h.01M14.5 15h.01" strokeWidth={3} />
        </>
      )}
      {theme === 'arithmetique' && <path d="M9 4l-2 16M17 4l-2 16M4 9h16M3.5 15h16" />}
      {theme === 'groupes' && (
        <>
          <path d="M18.5 9A7 7 0 0 0 6 7.5M5.5 15A7 7 0 0 0 18 16.5" />
          <path d="M6 4v3.5h3.5M18 20v-3.5h-3.5" />
        </>
      )}
      {theme === 'geometrie' && (
        <>
          <path d="M12 4l8 15H4z" />
          <circle cx="12" cy="13.6" r="3.4" />
        </>
      )}
      {theme === 'topologie-calcul-diff' && (
        <>
          <circle cx="12" cy="12" r="8" strokeDasharray="2.5 3" />
          <circle cx="12" cy="12" r="3.5" />
        </>
      )}
      {theme === 'probabilites' && (
        <>
          <rect x="4" y="4" width="16" height="16" rx="3.5" />
          <path d="M8.5 8.5h.01M15.5 15.5h.01M12 12h.01" strokeWidth={3} />
        </>
      )}
      {theme === 'methodes-numeriques' && <path d="M4 18h16M6 18l3-6 3 3 3-7 3 4" />}
    </svg>
  )
}

/** Chevron des volets (pointe vers le bas ; le CSS le retourne quand le volet est ouvert). */
export function Chevron({ taille = 18 }: { taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" {...TRAITS} strokeWidth={2}>
      <path d="M6 9l6 6 6-6" />
    </svg>
  )
}
