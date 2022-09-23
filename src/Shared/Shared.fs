namespace Shared

module Route =
    let builder typeName methodName = sprintf "/api/%s/%s" typeName methodName

type Depiction =
    | Filled
    | BallAndStick

type ViewBox = float * float * float * float

type Settings =
    { ViewBox: ViewBox option
      Depiction: Depiction
      ShowHydrogenAtoms: bool
      XRotation: float
      YRotation: float
      ZRotation: float }

type Assignment = { Settings: Settings; Sdf: string }

type ICinemolApi = { render: Assignment -> Async<string * string * ViewBox> }
