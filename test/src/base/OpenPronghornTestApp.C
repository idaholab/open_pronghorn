//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#include "OpenPronghornTestApp.h"
#include "OpenPronghornApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
OpenPronghornTestApp::validParams()
{
  InputParameters params = OpenPronghornApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

OpenPronghornTestApp::OpenPronghornTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  OpenPronghornTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

OpenPronghornTestApp::~OpenPronghornTestApp() {}

void
OpenPronghornTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  OpenPronghornApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"OpenPronghornTestApp"});
    Registry::registerActionsTo(af, {"OpenPronghornTestApp"});
  }
}

void
OpenPronghornTestApp::registerApps()
{
  registerApp(OpenPronghornApp);
  registerApp(OpenPronghornTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
OpenPronghornTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  OpenPronghornTestApp::registerAll(f, af, s);
}
extern "C" void
OpenPronghornTestApp__registerApps()
{
  OpenPronghornTestApp::registerApps();
}
