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

#pragma once

#include "MooseApp.h"

class OpenPronghornApp : public MooseApp
{
public:
  static InputParameters validParams();

  OpenPronghornApp(const InputParameters & parameters);
  virtual ~OpenPronghornApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);

  virtual std::string getInstallableInputs() const override;
};
