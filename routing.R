library(shiny.router)

source("uis/celloscope_ui.R")
# celloscope.page <- div(
#   titlePanel("celloscope"),
#   p("This is a contact page")
# )

source("uis/contact_ui.R")
# contact.page <- div(
#   titlePanel("Contact"),
#   p("This is a contact page")
# )

source("uis/legal_notice_ui.R")
# legal.notice.page <-  div(
#   titlePanel("Legal Notice"),
#   p("This is a default page")
# )

source("uis/privacy_policy_ui.R")
# privacy.policy.page <-  div(
#   titlePanel("Privacy Policy"),
#   p("This is a default page")
# )
  
  
router <- make_router(
  route("/", celloscope.page),
  route("contact", contact.page), 
  route("disclaimer", legal.notice.page),
  route("privacy", privacy.policy.page)
)

