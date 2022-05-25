privacy.policy.page <- shinyUI(
  fluidPage(
    p(
      includeHTML("www/privacy_policy.html")
    )
  )
)