// Temps et jours. Toute la logique reçoit l'instant présent en argument (horloge injectée) : jamais de Date.now() en dur
// hors de horlogeSysteme, pour que les tests tournent avec un faux temps.

/** Renvoie l'instant présent, en millisecondes depuis 1970. */
export type Horloge = () => number

export const horlogeSysteme: Horloge = () => Date.now()

const MS_JOUR = 86_400_000

/**
 * Jour civil local « AAAA-MM-JJ » de l'instant `t`. `decalageMinutes` = avance du fuseau sur UTC (+60 à Paris l'hiver) ;
 * par défaut celui de l'appareil à cet instant. Les échéances des cartes sont des jours : un changement de fuseau ne
 * déplace donc jamais une révision de plus d'un jour.
 */
export function jourLocal(t: number, decalageMinutes = -new Date(t).getTimezoneOffset()): string {
  return new Date(t + decalageMinutes * 60_000).toISOString().slice(0, 10)
}

const versMs = (jour: string): number => Date.parse(`${jour}T00:00:00Z`)

/** Le jour `n` jours après `jour` (n peut être négatif). */
export function ajouterJours(jour: string, n: number): string {
  return new Date(versMs(jour) + n * MS_JOUR).toISOString().slice(0, 10)
}

/** Nombre de jours de `a` à `b` (positif si `b` est après `a`). */
export function joursEntre(a: string, b: string): number {
  return Math.round((versMs(b) - versMs(a)) / MS_JOUR)
}
