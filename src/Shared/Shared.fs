namespace Shared


module Route =
    let builder typeName methodName = sprintf "/api/%s/%s" typeName methodName

type Settings =
    { ShowHydrogenAtoms: bool
      XRotation: float
      YRotation: float
      ZRotation: float }

type Assignment = { Settings: Settings; Sdf: string }

type ICinemolApi = { render: Assignment -> Async<string * string> }
