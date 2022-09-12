namespace Shared


module Route =
    let builder typeName methodName = sprintf "/api/%s/%s" typeName methodName

type Settings = { FilterHydrogens: bool }

type Assignment = { Settings: Settings; Sdf: string }

type ICinemolApi = { render: Assignment -> Async<string> }
