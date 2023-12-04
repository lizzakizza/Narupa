import os.path
import tempfile
import time
import warnings
import socket
import pytest
from cryptography import x509
from cryptography.x509.oid import NameOID
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.primitives import hashes
import datetime
from typing import Tuple
import grpc
from narupa.core import GrpcCredentials, GrpcServer
from narupa.trajectory import FrameServer, FrameClient, FrameData


def test_grpc_secure_connection_verification():
    """
    This unit test checks the validity and security of gRPC connections.

    It operates under the assumption that a connection is secure if:

     1. Data can be exchanged between an unsecured client and an unsecured server.
     2. Data can be exchanged between a secured client and a secured server.
     3. Data cannot be exchanged between a secured client and an unsecured server.
     4. Data cannot be exchanged between an unsecured client and a secured server.

    If all these conditions are met, the test passes, implying that the gRPC connections
    are secure and properly configured.

    Failures at different stages of the test can help pinpoint the nature of the issue:

     - Failure to establish an insecure connection may suggest a more general
       connection issue unrelated to SSL/TLS.
     - Failure to establish a secure connection could indicate a problem with
       the SSL/TLS implementation.
     - Failure to reject a partially secured connection (i.e., when only one
       party expects a secure connection) is critical and may suggest that
       encryption was never properly enforced in the first place.

    Note: Partially secured connections are expected to log an
    SSL_ERROR_SSL:WRONG_VERSION_NUMBER error.
    """

    def check(client_is_secure, server_is_secure):
        with ClientServer(client_is_secure, server_is_secure) as (client, server):
            return helper_client_server_communication(client, server)

    # Before checking if a secure connection can be established, first try to
    # establish an insecure connection. If the test fails here then there is
    # a more general connection issue at play here that is not necessary
    # caused directly by the SSL/TLS implementation.
    assert check(False, False), "Failed to establish an insecure connection"
    time.sleep(0.01)

    # Now check if a secured connection can be established. If the test fails
    # here it indicates that the error is associated specifically with the SSL/TLS
    # implementation.
    assert check(True, True), "Failed to establish a secure connection"
    time.sleep(0.01)

    # Now check if a partial connection can be established. This is a connection
    # in which either i) the client expects a secure connection but the server
    # does not, or ii) the server expects a secure connection the client does
    # not. This should always fail to form a connection. If the test fails here
    # then it indicates that encryption was never working in the first place.
    assert not check(True, False), "Critical: partial secure connection established"
    time.sleep(0.01)

    assert not check(False, True), "Critical: partial secure connection established"
    time.sleep(0.01)

    # Note that the last two checks do technically case an error and will thus
    # log an SSL_ERROR_SSL:WRONG_VERSION_NUMBER error.


def helper_client_server_communication(client, server):
    """Check that the specified server-client pair can communicate.

    This method takes a client-server pair and attempts to send data back and
    forth between them. If data can be sent and received then it is assumed
    that the connection is valid. If not it is assumed that the connection is
    invalid. The check is performed in this manner because errors associated
    with failed connection attempts are relegated to the background of Narupa
    at the moment. Thus, one cannot just test for exceptions upon attempting
    a connection.

    """
    trace = CallbackTrace()
    client.subscribe_frames_async(trace, 0.00)
    time.sleep(0.01)
    server.send_frame(0, simple_frame_data())
    time.sleep(0.01)
    return trace.success


def test_network_guard_functionality():
    """
    This unit test verifies the functionality and integrity of the network guard within :class:`GrpcServer`.

    Specifically, it checks the following scenarios:

     - Setting up an unsecured local server instance should not raise any
       issues or warnings.
     - Attempting to set up an unsecured non-local private server should raise
       a warning if the system is connected to a network.
     - Attempting to set up a public facing server without security should
       raise a warning.
     - No warnings should be issued if the server is secured, regardless of
       whether it is private or public.
     - No warnings should be issued if the `suppress_warnings` flag is set to
       True, regardless of the server's security status or network exposure.
    """

    local_ip = 'localhost'
    private_ip, connected_to_network = get_local_network_ip()
    open_ip = '[::]'
    port = 0

    with TemporaryCredentials(True) as (_, secure_credentials):
        insecure_credentials = GrpcCredentials(False)

        # Check that no issue is encountered when trying to set up an unsecure
        # local server instance.
        GrpcServer(address=local_ip, port=port, credentials=insecure_credentials).close()

        # Check that a warning is issued when setting up a private but non-
        # local server. Not that this is, and only can be, performed if the
        # system is connected to a network.
        if connected_to_network:
            with pytest.warns(UserWarning, match=r"(?<=Security Warning:).+"):
                GrpcServer(address=private_ip, port=port, credentials=insecure_credentials).close()

        # Check that a warning is issued when possibly creating a public facing
        # server without security.
        with pytest.warns(UserWarning, match=r"(?<=Security Warning:).+"):
            GrpcServer(address=open_ip, port=port, credentials=insecure_credentials).close()

        # Check that the warnings are not issued when the server is secured.
        with warnings.catch_warnings():
            warnings.filterwarnings("error", message=r"Security Warning:")
            if connected_to_network:
                GrpcServer(address=private_ip, port=port, credentials=secure_credentials).close()

            GrpcServer(address=open_ip, port=port, credentials=secure_credentials).close()

        # Check that the warnings are not issued when the suppress_warnings is
        # set.
        with warnings.catch_warnings():
            warnings.filterwarnings("error", message=r"Security Warning:")
            if connected_to_network:
                GrpcServer(address=private_ip, port=port,
                           credentials=insecure_credentials,
                           suppress_warnings=True).close()

            GrpcServer(address=open_ip, port=port,
                       credentials=insecure_credentials,
                       suppress_warnings=True).close()


def test_grpc_credentials_integrity():
    """
    This unit test verifies the integrity and functionality of the :class:`GrpcCredentials` class.

    It covers the following cases:

     - Creation of secure and insecure credentials.
     - Correct behaviour of the equality operator for secure and insecure credential objects.
     - Serialization and deserialization of both secure and insecure credentials, ensuring they
       remain the unchanged following this process.
     - Correct identification and handling of invalid or inappropriate parameters when
       validating the credentials, such as missing files or mismatched attributes.
     - Assertion of proper object types when invoking methods to retrieve server and client
       credentials from a secure :class:`GrpcCredentials` object.
     - Verification that an error is correctly raised when trying to retrieve server and client
       credentials from an insecure :class:`GrpcCredentials` object.
    """

    # Check that secure and insecure credentials can be created without issue
    with TemporaryCredentials(True) as (secure_credentials_a, secure_credentials_b):
        insecure_credential = GrpcCredentials(False)

        # Check that the equality operator behaves as intended
        assert insecure_credential == insecure_credential,   "Equality check 1"
        assert secure_credentials_a == secure_credentials_a, "Equality check 2"
        assert secure_credentials_a != insecure_credential,  "Equality check 3"
        assert secure_credentials_a != secure_credentials_b, "Equality check 4"

        # Check that credentials can be serialised and deserialised without
        # mutation.

        # Check that the insecure credential can be serialised and deserialised.
        with tempfile.TemporaryDirectory() as temp_dir:
            path_u = os.path.join(temp_dir, 'temp_u.json')
            insecure_credential.dump_grpc_credentials_to_json_file(path_u)
            insecure_credential_reloaded = GrpcCredentials.load_grpc_credentials_from_json_file(path_u)
            assert insecure_credential == insecure_credential_reloaded, "Failed to (de)serialised unsecure credentials"

            path_s = os.path.join(temp_dir, 'temp_s.json')
            secure_credentials_a.dump_grpc_credentials_to_json_file(path_s)
            secure_credentials_a_reloaded = GrpcCredentials.load_grpc_credentials_from_json_file(path_s)

            assert secure_credentials_a == secure_credentials_a_reloaded, "Failed to (de)serialised secure credentials"

        # Check that the validation methods catches missing files or mismatching attributes

        # If the credentials are secure then file paths must be specified
        with pytest.raises(ValueError):
            insecure_credential.secure = True
            insecure_credential.validate()

        # If file paths are given but do not point to valid files
        with pytest.raises(FileNotFoundError):
            secure_credentials_a.certificate_file = "an/invalid/path.json"
            secure_credentials_a.validate()

        # If credentials are insecure then it makes no sense to provide file paths
        with pytest.raises(ValueError):
            secure_credentials_a.secure = False
            secure_credentials_a.validate()

        # Check that server and client credentials can be returned by the
        # "get" methods. There is no need to test for the validity of the
        # returns as this will show up in tests elsewhere. But first make
        # sure that an error is thrown if the user attempts to perform this
        # action on an unsecure credential set.
        with pytest.raises(AssertionError):
            GrpcCredentials(False).get_client_credentials()

        with pytest.raises(AssertionError):
            GrpcCredentials(False).get_server_credentials()

        assert isinstance(secure_credentials_b.get_client_credentials(), grpc.ChannelCredentials)
        assert isinstance(secure_credentials_b.get_server_credentials(), grpc.ServerCredentials)


class TemporaryCredentials:
    """Context manager to create temporary certificate files in a dedicated directory.

    When exiting the context, the directory and all its contents are
    automatically cleaned up.

    The security credentials generated through the use of this class are to be
    used for testing only! They should never be used in production.

    :param secure: indicates if credentials should be secure
    :param client_cn: common name of the client
    :param server_cn: common name of the server

    :Example:

    with TemporaryCredentials() as (client_credentials, server_credentials):
        # Use the GrpcCredentials instances credentials
        ...

    # Here, the directory and all its files have been deleted
    """

    def __init__(self, secure: bool = True, client_cn: str = 'client',
                 server_cn: str = 'localhost'):

        self.secure = secure
        self.client_cn = client_cn
        self.server_cn = server_cn

        self.temp_dir = None

    def __enter__(self):
        """Called when entering the context (i.e., the `with` block).

        :return: credential objects for the client and server.
        :rtype: Tuple[GrpcCredentials, GrpcCredentials]
        """
        if self.secure:
            self.temp_dir = tempfile.TemporaryDirectory()
            self.path = self.temp_dir.name
            return self._generate_certificates(self.client_cn, self.server_cn)
        else:
            return GrpcCredentials(False), GrpcCredentials(False)

    def __exit__(self, *args):
        """Called when exiting the context (i.e., the `with` block).
        """
        if self.temp_dir is not None:
            self.temp_dir.cleanup()

    def _generate_certificates(self, client_cn, server_cn):
        """Generate the CA, server, and client private keys and certificates.

        This method calls the helper methods `_generate_ca`,
        `_generate_key_and_csr`, and `_sign_csr` to generate the necessary
        keys and certificates. The generated files are saved in the temporary
        directory.


        :return: credential objects for the client and server.
        :rtype: Tuple[GrpcCredentials, GrpcCredentials]
        """

        ca_key, ca_cert = self._generate_ca()
        client_key, client_csr = self._generate_key_and_csr(client_cn)
        server_key, server_csr = self._generate_key_and_csr(server_cn)

        # sign CSRs with the CA key
        self._sign_csr(client_csr, ca_key, ca_cert)
        self._sign_csr(server_csr, ca_key, ca_cert)

        # Create and return the client and server `GrpcCredentials` instances
        def jp(target):
            return os.path.join(self.path, target)

        client_creds = GrpcCredentials(
            True,
            jp("ca.crt"),
            jp(f"{client_cn}.crt"),
            jp(f"{client_cn}.key"))

        server_creds = GrpcCredentials(
            True,
            jp("ca.crt"),
            jp(f"{server_cn}.crt"),
            jp(f"{server_cn}.key"))

        return client_creds, server_creds

    def _generate_key_and_csr(self, common_name):
        """Generate a private key and a Certificate Signing Request (CSR).

        :param common_name: The common name to use in the CSR
        :type common_name: str
        """
        key = rsa.generate_private_key(
            public_exponent=65537,
            key_size=2048,
        )

        csr = x509.CertificateSigningRequestBuilder().subject_name(x509.Name([
            x509.NameAttribute(NameOID.COMMON_NAME, common_name),
        ])).sign(key, hashes.SHA256())

        with open(os.path.join(self.path, f"{common_name}.key"), "wb") as f:
            f.write(key.private_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PrivateFormat.TraditionalOpenSSL,
                encryption_algorithm=serialization.NoEncryption(),
            ))

        with open(os.path.join(self.path, f"{common_name}.csr"), "wb") as f:
            f.write(csr.public_bytes(serialization.Encoding.PEM))

        return key, csr

    def _generate_ca(self):
        """Generate a Certificate Authority (CA), private key and self-signed certificate.
        """
        ca_key = rsa.generate_private_key(
            public_exponent=65537,
            key_size=2048,
        )

        ca_cert = (
            x509.CertificateBuilder()
            .subject_name(x509.Name([x509.NameAttribute(NameOID.COMMON_NAME, u"")]))
            .issuer_name(x509.Name([x509.NameAttribute(NameOID.COMMON_NAME, u"")]))
            .public_key(ca_key.public_key())
            .serial_number(x509.random_serial_number())
            .not_valid_before(datetime.datetime.utcnow())
            .not_valid_after(datetime.datetime.utcnow() + datetime.timedelta(days=365))
            .add_extension(x509.BasicConstraints(ca=True, path_length=None), critical=True)
            .sign(ca_key, hashes.SHA256())
        )

        with open(os.path.join(self.path, "ca.crt"), "wb") as f:
            f.write(ca_cert.public_bytes(serialization.Encoding.PEM))

        with open(os.path.join(self.path, "ca.key"), "wb") as f:
            f.write(ca_key.private_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PrivateFormat.TraditionalOpenSSL,
                encryption_algorithm=serialization.NoEncryption(),
            ))

        return ca_key, ca_cert

    def _sign_csr(self, csr, ca_key, ca_cert):
        """Sign a Certificate Signing Request (CSR) with the Certificate Authority (CA) key.

        :param csr: The Certificate Signing Request to be signed
        :param ca_key: The private key of the Certificate Authority
        :param ca_cert: The certificate of the Certificate Authority
        """
        cert = (
            x509.CertificateBuilder()
            .subject_name(csr.subject)
            .issuer_name(ca_cert.subject)
            .public_key(csr.public_key())
            .serial_number(x509.random_serial_number())
            .not_valid_before(datetime.datetime.utcnow())
            .not_valid_after(datetime.datetime.utcnow() + datetime.timedelta(days=365))
            .add_extension(x509.BasicConstraints(ca=False, path_length=None), critical=True)
            .sign(ca_key, hashes.SHA256())
        )
        with open(os.path.join(
                self.path,
                f"{cert.subject.get_attributes_for_oid(NameOID.COMMON_NAME)[0].value}.crt"),
                "wb") as f:
            f.write(cert.public_bytes(serialization.Encoding.PEM))


def simple_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.arrays["indices"] = [0, 1, 3]
    basic_frame_data.values["string"] = "str"
    basic_frame_data.values["bool"] = False
    return basic_frame_data


class ClientServer:
    """Use to quickly produce a client server pair"""

    def __init__(self, secure_client, secure_server):
        self._secure_client, self._secure_server = secure_client, secure_server
        self._client, self._server = None, None

    def __enter__(self):
        # Generate secure credentials
        with TemporaryCredentials() as (client_credentials, server_credentials):
            # Overwrite with insecure credentials as needed
            server_credentials = server_credentials if self._secure_server else GrpcCredentials(False)
            client_credentials = client_credentials if self._secure_client else GrpcCredentials(False)
            # Start up the client and server instances
            self._server = FrameServer(
                address='localhost', port=0, credentials=server_credentials)
            self._client = FrameClient.establish_channel(
                address='localhost', port=self._server.port, credentials=client_credentials)
            # Return the client and server
            return self._client, self._server

    def __exit__(self, *args):
        # Upon exiting the context manager the client and server are shut down
        self._client.close(), self._server.close()


class CallbackTrace:
    """Used to check if data is transmitted"""
    success = False
    def __call__(self, *args, **kwargs): self.success = True


def get_local_network_ip():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        # doesn't have to be reachable
        s.connect(('10.255.255.255', 1))
        ip = s.getsockname()[0]
        connected_to_network = True
    except:
        ip = None
        connected_to_network = False
    finally:
        s.close()
    return ip, connected_to_network
