from typing import Optional, Union
from os.path import isfile
import json
import grpc


class GrpcCredentials:
    """Manages the credentials necessary for gRPC communication.

    This class is used to determine whether a user wishes to establish a secure
    and encrypted connection, and offers functionally to facilitate this. If
    using encryption, this entity will store paths to the required certifications
    and keys. The `get_server_credentials` and `get_client_credentials` methods
    can then be invoked to get the necessary `grpc.ChannelCredentials` entity.
    For insecure connections, the `secure` attribute is set to `False` which
    is indicates to the server/client that the user expressly wishes to use
    an unencrypted and insecure channel.

    Note that certificate and key files are only required for secure connections;
    i.e. when `secure=True`.

    :param secure: a boolean indicating whether a gRPC connection should be secured.
    :param root_certificate_authority_file: path to the root certificate authority file.
    :param certificate_file: path to the client or server's certificate file.
    :param key_file: path to the client or server's private key file.

    """
    def __init__(self, secure: bool,
                 root_certificate_authority_file: Optional[str] = None,
                 certificate_file: Optional[str] = None,
                 key_file: Optional[str] = None):

        self.secure = secure
        self.root_certificate_authority_file = root_certificate_authority_file
        self.certificate_file = certificate_file
        self.key_file = key_file

        self.validate()

    def get_server_credentials(self) -> grpc.ChannelCredentials:
        """Construct secure server credentials.

        :return: a `grpc.ChannelCredentials` instance that can be used by a gRPC
            server for securing and authenticating connections.
        """

        # Calling `get_server_credentials` on an insecure `GrpcCredentials`
        # instance is undefined.
        assert self.secure, "get_server_credentials is only valid for secure servers"

        # Load in the certificate and key files
        root_certificate_authority, certificate, key = self.__load_files()

        # Parse them into a `grpc.ChannelCredentials` instance and return.
        return grpc.ssl_server_credentials(
            [(key, certificate)],
            root_certificates=root_certificate_authority,
            require_client_auth=True
        )

    def get_client_credentials(self) -> grpc.ChannelCredentials:
        """Construct secure client credentials.

        :return: a `grpc.ChannelCredentials` instance that can be used by a gRPC
            client for securing and authenticating connections.
        """

        # Calling `get_client_credentials` on an insecure `GrpcCredentials`
        # instance is undefined.
        assert self.secure, "get_client_credentials is only valid for secure clients"

        # Load in the certificate and key files
        root_certificate_authority, certificate, key = self.__load_files()

        # Parse them into a `grpc.ChannelCredentials` instance and return.
        return grpc.ssl_channel_credentials(
            private_key=key, certificate_chain=certificate,
            root_certificates=root_certificate_authority
        )

    def validate(self):
        """Validate credential setup.

        This method ensures consistency between the desired level of security,
        as specified by `secure`, and the files provided. Additional checks
        are also performed to check that all file paths are valid.

        :raises:
            ValueError: thrown for secure connections when the credentials are invalid or missing.
            FileNotFoundError: thrown for secure when a specified file does not exist.
            ValueError: thrown when credentials are provided for non-secure connections.
        """

        if self.secure:
            # If secure credentials are to be used then check that paths for
            # all required keys and certificates have been specified
            if not all(map(_is_valid_path, (
                    self.root_certificate_authority_file,
                    self.certificate_file, self.key_file))):
                raise ValueError("For secure connections, certificate and key files must not be null.")

            # and that the file paths are valid.
            if not all(map(isfile, (
                    self.root_certificate_authority_file,
                    self.certificate_file, self.key_file))):
                raise FileNotFoundError("One or more certificate or key file paths are invalid.")

            # It might be worth adding certificate validation later on down
            # the line, however, the current tests will suffice for now.
        else:
            # For insecure credentials, no files are required; an exception is
            # still thrown as the user may have intended to use secure credentials.
            if any(map(_is_valid_path, (
                    self.root_certificate_authority_file,
                    self.certificate_file, self.key_file))):
                raise ValueError("Certificate and key files are provided but secure connection is disabled.")

    def __load_files(self):
        """Read in and return the certificate and key files."""
        with open(self.root_certificate_authority_file, 'rb') as file_handle:
            root_certificate_authority = file_handle.read()

        with open(self.certificate_file, 'rb') as file_handle:
            certificate = file_handle.read()

        with open(self.key_file, 'rb') as file_handle:
            key = file_handle.read()

        return root_certificate_authority, certificate, key

    def __eq__(self, other: 'GrpcCredentials') -> bool:
        # Only equality checks against other `GrpcCredentials` instance are
        # supported.
        assert isinstance(other, self.__class__)
        self.validate()
        other.validate()
        if self.secure and other.secure:
            # If both credentials are secure then check if the contents of
            # their key and certificate files match.
            return all([i == j for i, j in zip(
                self.__load_files(), other.__load_files())])
        elif not self.secure and not other.secure:
            # If both are insecure then they will match, assuming that they
            # are both valid.
            return True
        else:
            # If one is secure and the other is not, then they do not match.
            return False

    def __ne__(self, other: 'GrpcCredentials') -> bool:
        assert isinstance(other, self.__class__)
        return not self == other

    @classmethod
    def load_grpc_credentials_from_json_file(cls, path) -> 'GrpcCredentials':
        """Instantiate a `GrpcCredentials` instance from a JSON structured file.

        :param path: path to a JSON structured file storing a serialised
            `GrpcCredentials` instance.

        """
        with open(path) as file_handle:
            data = json.load(file_handle)
            return cls(
                data["Secure"], data["RootCertificateAuthorityFile"],
                data["CertificateFile"], data["KeyFile"])

    def dump_grpc_credentials_to_json_file(self, path):
        """Write out a `GrpcCredentials` instance into a JSON structured file.

        :param path: path to the JSON structured file in which the serialised
            `GrpcCredentials` instance is to be stored.

        """
        self.validate()

        with open(path, 'w') as file_handle:
            file_handle.write(
                json.dumps(
                    {
                        "Secure": self.secure,
                        "RootCertificateAuthorityFile": self.root_certificate_authority_file,
                        "CertificateFile": self.certificate_file,
                        "KeyFile": self.key_file
                    }
                    , indent=4))


def _is_valid_path(path: Union[str, None]) -> bool:
    """Tests if a string is a valid non-blank, non-whitespace string.

    This is used to check whether a string is actually intended to represent a
    path. Path attributes of `GrpcCredentials` instances will commonly be set
    to `None` if they are neither specified nor required. However, credentials
    are stored in a JSON file the paths will commonly be set to an empty string,
    i.e `""`. This is due to the inability of the complimentary C-sharp
    `GrpcCredentials` class's inability to store nullable strings. As it is
    currently based on Unity 2019.X which uses a deprecated version C-sharp
    compiler. This can safely be removed once Narupa is migrated to a modern
    version of Unity.

    :param path: possibly empty path to a key or certificate file.

    :return: A boolean indicating that ``path`` is a valid path string.


    """
    # Check `path` is 1) a string, 2) not empty, & 3) not whitespace
    return isinstance(path, str) and (path and not path.isspace())
