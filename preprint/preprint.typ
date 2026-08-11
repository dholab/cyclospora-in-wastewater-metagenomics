#import "@preview/arkheion:0.1.2": arkheion

#set terms(hanging-indent: 1.5em)
#set table(inset: 5pt, stroke: none)
#show table: set text(size: 8.5pt)
#show table: set par.line(numbering: none)
#show link: underline

#let review-lines(body) = {
  set par.line(
    numbering: number => text(size: 7pt, fill: luma(45%))[#number],
    number-align: right,
    number-margin: left,
    number-clearance: 6pt,
    numbering-scope: "document",
  )
  body
}

#let horizontalrule = line(start: (25%, 0%), end: (75%, 0%))
#let divider = if "divider" in std { divider } else { horizontalrule }

#show figure.where(kind: table): set figure.caption(position: top)
#show figure.where(kind: image): set figure.caption(position: bottom)
#show figure.caption: caption => align(left, caption)

#let manuscript-authors = (
$for(authors)$
  (name: "$authors.name$"),
$endfor$
)

#let corresponding-email = "$corresponding.email_user$" + "@" + "$corresponding.email_domain$"

#show: arkheion.with(
  title: "$title$",
  authors: manuscript-authors,
  custom-authors: [
    #align(center)[
      $for(authors)$*$authors.name$*#super[$authors.marks$]$sep$, $endfor$
      #v(0.6em)
      #set text(size: 9pt)
      $for(affiliations)$#super[$affiliations.mark$] $affiliations.text$ \
      $endfor$
      #super[\*] Corresponding author: $corresponding.name$ · #link("mailto:" + corresponding-email)[#corresponding-email]
    ]
  ],
  abstract: review-lines([
$abstract$
  ]),
)

#review-lines[
$body$
]
