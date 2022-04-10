## workHORSE is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - workHORSE: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## workHORSE is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

tabPanel(
  "Welcome",
  h1(
    "workHORSE: working Health Outcomes Research Simulation Environment"
  ),
  br(),
  br(),

  p(
    "Welcome to the prototype workHORSE web-based decision support tool designed to model the potential for equitable health gain and cost-effectiveness of the NHS Health Check Programme.",
    style = "font-size:16px"
  ),
  br(),
  p(
    "At the core of workHORSE there is a microsimulation model that tracks the health outcomes of synthetic individuals under different implementation scenarios of NHS Health Check Programme, locally and nationally. The synthetic individuals are coming from a validated 'close-to-reality' synthetic population primed from the Health Survey for England series that closely resembles the health behaviours and exposures of the English population.",
    style = "font-size:16px"
  ),
  br(),
  p(
    "The user can design different NHS Health Check Programme implementation scenarios and evaluate and compare their effectiveness, equity, and cost-effectiveness. workHORSE generates rich outputs that could support real world decision making regarding the  NHS Health Check Programme. For best results, the user should inform the scenarios builder with local data regarding the existing implentation of NHS Health Check Programme in their area of interest.",
    style = "font-size:16px"
  ),
  br(),
  p(
    "This tool has been co-produced with a wide group of stakeholders with expertise in the commissioning and provision of the NHS Health Checks Programme, third sector organisations and national organisations. We are grateful for all their invaluable feedback and creative ideas that made a real impact on the usability and user-friendlyness of workHORSE.",
    style = "font-size:16px"
  ),
  br(),
  p(
    "WorkHORSE was a 2-year project funded by the National Institute of Health Research (NIHR) NIHR health technology assessment programme (project HTA 16/165/01). The funders had no role in study design, data collection and analysis, decision to publish, or development of workHORSE.",
    style = "font-size:16px"
  ),
  br(),
  p(
    "For more information about workHORSE please contact Chris Kypridemos (ckyprid  at liverpool.ac.uk) or visit: ",
    style = "font-size:16px;display:inline"),
    a(href = "https://www.workhorsesim.net/", "www.workHORSEsim.net",
      style = "font-size:16px", target = "_blank"),
  icon = icon("door-open")
)
