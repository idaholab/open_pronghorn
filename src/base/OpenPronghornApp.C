#include "OpenPronghornApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
OpenPronghornApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

OpenPronghornApp::OpenPronghornApp(InputParameters parameters) : MooseApp(parameters)
{
  OpenPronghornApp::registerAll(_factory, _action_factory, _syntax);
}

OpenPronghornApp::~OpenPronghornApp() {}

void
OpenPronghornApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<OpenPronghornApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"OpenPronghornApp"});
  Registry::registerActionsTo(af, {"OpenPronghornApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
OpenPronghornApp::registerApps()
{
  registerApp(OpenPronghornApp);
  ModulesApp::registerApps();
}

std::string
OpenPronghornApp::getInstallableInputs() const
{
  return OPEN_PRONGHORN_INSTALLABLE_DIRS;
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
OpenPronghornApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  OpenPronghornApp::registerAll(f, af, s);
}
extern "C" void
OpenPronghornApp__registerApps()
{
  OpenPronghornApp::registerApps();
}
