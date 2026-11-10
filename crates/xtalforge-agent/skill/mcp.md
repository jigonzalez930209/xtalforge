# MCP client setup

`xtalforge mcp` is a Model Context Protocol server over stdio. It speaks
JSON-RPC 2.0, newline-delimited, and implements protocol version `2024-11-05`.

## Client configuration

A client that starts the command itself uses a configuration like:

```json
{
  "mcpServers": {
    "xtalforge": {
      "command": "xtalforge",
      "args": ["mcp"]
    }
  }
}
```

No environment variables are required. The server uses its standard input and
output only and writes no logs there.

## Tools

`tools/list` returns one tool per verb. The name is the verb
(`load`, `info`, `supercell`, ...) and `inputSchema` is generated from the same
catalogue the CLI and bindings use. Required arguments are listed in
`required`; defaults are in the schema. Unknown arguments are rejected.

## Results

A successful call returns:

```json
{
  "content": [{"type": "text", "text": "{\"kind\": \"info\", ...}"}],
  "isError": false
}
```

A failed call returns `isError: true` and the text payload is
`{"code": "...", "message": "..."}`. The server never terminates on a bad tool
call.

## Lifecycle

- `initialize` returns `protocolVersion`, `capabilities.tools`,
  `serverInfo` and `instructions`.
- `notifications/initialized` is accepted without a response.
- `ping` returns `{}`.
- `notifications/cancelled` is accepted. Verbs are synchronous, so cancellation
  is applied between calls; to stop a long verb, terminate the process.
- An unknown method returns JSON-RPC error `-32601`; malformed JSON returns
  `-32700`.

## Permissions

The server can read and write any file the process can. Run it under the user
account that should own the documents; there is no sandbox beyond the size
limits.
