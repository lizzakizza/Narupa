# Secure Client-Server Communication
This guide provides an overview of and a walkthrough for establishing secure client-server communication within Narupa.
Secure communication is essential to maintain data privacy and prevent unauthorized access, particularly when dealing with sensitive data.

Narupa leverages gRPC, a modern open-source high-performance Remote Procedure Call (RPC) framework, to handle communication between client and server.
However, it's important to note that such communication is **not secure or encrypted by default**.

This lack of security is commonly not an issue when running both the client and server on the same machine, in a local isolated environment.
However, this can present as significant risk if the client and server are not coterminous or are exposed to external network traffic.
Without appropriate safeguards, data can be intercepted, manipulated, or stolen.
Therefore, Narupa blocks attempts to set up a server that is open to external network traffic due to the inherent security risks.

To combat these risks, when dealing with remote servers, both the client and the server must undergo an authentication process.
Following, which all data transmitted is encrypted.

This guide walks through the steps required to i) generate a Certificate Authority (CA), and the necessary keys and certificates for client and server, and ii) pass the newly generated credentials to the client and server
This process ensures secure SSL/TLS encrypted communication, preventing any potential eavesdropping or data tampering.

## Generating Keys for Secure Communication
This section provides the steps to generate a Certificate Authority (CA), as well as the keys and certificates required for a client and a server to establish secure and encrypted SSL/TLS communication.

### Certificate Authority
Start by generating a private key, `ca.key`, for the CA:
```bash
openssl genpkey -algorithm RSA -out ca.key
```
Then create a self-signed certificate, `ca.crt`, for the CA.
You may leave the requested certificate information blank.
```bash
openssl req -new -x509 -days 365 -key ca.key -out ca.crt
```

### Client Key and Certificate
Generate a private key, `client.key`, for the client:
```bash
openssl genpkey -algorithm RSA -out client.key
```
Create a certificate signing request, `client.csr`, for the client. Most of the information prompts may be left blank.
However, a Common Name (CN) must be provided.
The choice of name is somewhat arbitrary, but leaving this field blank will result in an invalid certificate.
```bash
openssl req -new -key client.key -out client.csr
```
Sign the client certificate to generate the certificate file `client.crt`:
```bash
openssl x509 -req -in client.csr -CA ca.crt -CAkey ca.key -CAcreateserial -out client.crt -days 365
```

### Server Key and Certificate
Generate a private key, `server.key`, for the server:
```bash
openssl genpkey -algorithm RSA -out server.key
```
Create a certificate signing request, `server.csr`, for the server.
It's important to provide a correct **Common Name (CN)** when prompted:
```bash
openssl req -new -key server.key -out server.csr
```
The Common Name (CN) is an **important** part of a certificate used in secure gRPC network connections.
It should match the server's hostname, specifically the name that clients will be using to connect.
The client checks the CN against the server's address. If they don't match, the client rejects the connection.
This protects against "man-in-the-middle" attacks.
You can use wildcards in the CN (e.g. `*.example.com` will accept any subdomain of `example.com`).
If you're working locally, set the CN to `localhost` as this is the default broadcast address for Narupa servers.
Finally, sign the server certificate:
```bash
openssl x509 -req -in server.csr -CA ca.crt -CAkey ca.key -CAcreateserial -out server.crt -days 365
```
Note that the common name **cannot** be an IP address, as this is strictly [prohibited](https://github.com/grpc/grpc/blob/7910554cdd9ab4fa1eae77aa16739918ba572fea/src/core/tsi/ssl_transport_security.cc#L2326) by gRPC.
However, if one must connect to a server via an IP address, rather than a fully qualified domain name, then the "*server's name*" can be set to an IP address using a Subject Alternative Name (SAN).
This is done by creating a `san.ext` file like so:
```text
subjectAltName = @alt_names

[alt_names]
DNS.1 = myserver.com
IP.1 = 192.168.1.1
```
and then passing the `san.ext` file into the server-side certificate generation by appending `-extfile san.ext` to the command.

### Notes
The validity of the client and server certificates may be checked against the root certificate-authority by invoking the following commands:
```bash
openssl verify -CAfile ca.crt client.crt
openssl verify -CAfile ca.crt server.crt
```

## Storage and Deployment of Keys
It is vital to ensure that the keys and certificates are kept secure and private once they have been generated.

### Storage
The client and server require access to their respective key/certificate pairs and the CA's certificate, specifically:
- The client machine needs `client.key`, `client.crt`, and `ca.crt`.
- The server machine requires `server.key`, `server.crt`, and `ca.crt`.

To form a secure connection, the client and server must each be provided with an appropriate set of credentials via a configuration file.
This JSON-structured file indicates whether communication should be secure and where to locate the necessary key and certificate files.

A complete configuration file includes four key-value pairs:
- `Secure`: Boolean indicating if a secure connection is necessary.
- `RootCertificateAuthorityFile`: Absolute path to the root CA's certificate file (i.e., `ca.crt`).
- `CertificateFile`: Absolute path to the client or server's certificate file (`client.crt` or `server.crt`).
- `KeyFile`: Absolute path to the client or server's key file (`client.key` or `server.key`).


While a credentials configuration file is only needed for a secure connection, you can still use it with `Secure` set to `false` and file paths set as empty strings.
Although this makes the configuration file somewhat redundant, it is supported as the JSON file serves as a serialized C#/Python `GrpcCredentials` instance.

An example credentials configuration file has been provided below:
```json
{
    "Secure": true,
    "RootCertificateAuthorityFile": "PATH/TO/ca.crt",
    "CertificateFile": "PATH/TO/server.crt",
    "KeyFile": "PATH/TO/server.key"
}
```
### Deployment

#### Server
For Python based Narupa servers, the credentials are parsed into a `GrpcCredentials` instance which is then passed by keyword into a supported runner, like so:
```Python
from narupa.openmm import serializer
from narupa.ase.openmm import ASEOpenMMRunner
from narupa.core import GrpcCredentials

# Read the neuraminidase xml file, which defines the simulation parameters
with open('PATH/TO/neuraminidase.xml') as infile:
    simulation = serializer.deserialize_simulation(infile.read())

# Load in the credentials configuration file
credentials = GrpcCredentials.load_grpc_credentials_from_json_file(
    'PATH/TO/server_credentials.json')

# Create the server and provide it with the required connection credentials.
server = ASEOpenMMRunner(simulation, credentials=credentials)

# Finally run the server
server.run()
```

#### Client
The credential selection process for clients, such as [Narupa-iMD](https://gitlab.com/intangiblerealities/narupa-applications/narupa-imd) , is still a work in progress and is thus subject to change.
Currently, during startup, the clients scan Unity's persistent data directory (`Application.persistentDataPath`) for a `defaults.json` file in the `Credentials` subdirectory.
If a valid configuration file is found, it is used to set up a secure connection.

The exact location of Unity's persistent data directory varies based on certain factors, as outlined [here](https://docs.unity3d.com/ScriptReference/Application-persistentDataPath.html).
However, the typical paths for the three most common operating systems are provided below:

- **Windows Editor and Standalone Player**: `%userprofile%\AppData\LocalLow\<companyname>\<productname>`.
- **Linux**: `$XDG_CONFIG_HOME/unity3d` or `$HOME/.config/unity3d`.
- **Mac**: commonly points to the user's Library folder, often hidden. In recent Unity releases, user data is written into `~/Library/Application Support/company name/product name`. Older Unity versions wrote into the `~/Library/Caches` folder or `~/Library/Application Support/unity.company name.product name`.

These directories are automatically created upon running the Narupa client for the first time.
The placeholders `<companyname>` and `<productname>` will resolve to `Intangible Realities Laboratory` and `Narupa iMD` respectively for the Narupa-iMD client.
Therefore, for a Windows system, the credentials configuration file would be located at: `%userprofile%\AppData\LocalLow\Intangible Realities Laboratory\Narupa iMD\Credentials\defaults.json`.