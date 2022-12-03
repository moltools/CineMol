[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Maintainer](https://img.shields.io/badge/Maintainer-davidmeijer-blue)](https://github.com/davidmeijer)
[![Generic badge](https://img.shields.io/badge/Version-alpha-green.svg)](https://shields.io/)

# CineMol

<img src="https://github.com/davidmeijer/cinemol/blob/main/logo.png" alt="logo" width="100">

A [web application](https://moltools.nl/cinemol) for drawing SVG images of small molecules!

:warning: This project is a work in progress :warning:

## Install pre-requisites

You'll need to install the following pre-requisites in order to build CineMol:

* [.NET Core SDK](https://www.microsoft.com/net/download) 6.0 or higher
* [Node 16](https://nodejs.org/en/download/)

## Starting the application

Before you run the project **for the first time only** you must install dotnet "local tools" with this command:

```bash
dotnet tool restore
```

To run the application in watch mode use the following command:

```bash
dotnet run
```

Then open `http://localhost:8080` in your browser.

The build project in root directory contains a couple of different build targets. You can specify them after `--` (target name is case-insensitive).

To run client tests in watch mode (you can run this command in parallel to the previous one in new terminal):

```bash
dotnet run -- RunTests
```

Client tests are available under `http://localhost:8081` in your browser.

## Bundle the application

Create production bundle:

```bash
dotnet run Bundle
```

You can find the production bundle in the newly created ```deploy``` folder.
