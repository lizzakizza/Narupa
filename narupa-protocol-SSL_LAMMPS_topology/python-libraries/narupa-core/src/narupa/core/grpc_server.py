# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing a wrapper around the running of GRPC servers.
"""
import logging
import textwrap
import warnings
from concurrent import futures
import ipaddress
import socket

from typing import Optional, Tuple

import grpc

from narupa.core.grpc_credentials import GrpcCredentials

DEFAULT_SERVE_ADDRESS = 'localhost'
DEFAULT_CONNECT_ADDRESS = 'localhost'

# We expect that reserving a large number of threads should not present a
# performance issue. Each concurrent GRPC request requires a worker, and streams
# occupy those workers indefinitely, so several workers must be available for
# each expected client.
DEFAULT_MAX_WORKERS = 1000


class GrpcServer:
    """
    A base class for running GRPC servers that handles the starting and closing
    of the underlying server.

    :param address: The IP address at which to run the server.
    :param port: The port on which to run the server.
    :param credentials: Credentials specifying whether the server should be secured.
    """

    def __init__(
            self,
            *,
            address: str,
            port: int,
            credentials: Optional[GrpcCredentials] = None,
            **kwargs):

        # Future iterations of GrpcServer will require credentials to be
        # provided.
        if credentials is None:
            # In later versions a deprecation warning will be issued. But for
            # now the credentials will just be default initialised.
            # warnings.warn(
            #     "This call method is now deprecated please ensure that GrpcServer"
            #     " instances are provided with credentials during instantiation.",
            #     DeprecationWarning)

            credentials = GrpcCredentials(False)

        # Before potentially opening up an unsafe port perform some safety checks.
        _network_guard(
            address, port, credentials.secure,
            kwargs.get('suppress_warnings', False))

        # Set up the server and associated services
        self.server = grpc.server(
            futures.ThreadPoolExecutor(max_workers=kwargs.get("max_workers", DEFAULT_MAX_WORKERS)),
            options=(
                # Prevent multiple servers from hosting on the same port
                ('grpc.so_reuseport', 0),
            ))

        self.setup_services()

        # Open up a port for accepting RPCs. Check the `credentials` entity to
        # see whether the port should be encrypted and secured. On occasion,
        # users may set the port number to zero. This instructs gRPC to use the
        # first available port. The "add port" methods return the port number,
        # hence the `_port` attribute is assigned in this manner.
        if credentials.secure:
            self._port = self.server.add_secure_port(
                f"{address}:{port}", credentials.get_server_credentials())
        else:
            self._port = self.server.add_insecure_port(
                address=f"{address}:{port}")

        self._address = address


        # Enable logging
        self.logger = logging.getLogger(__name__)
        self.logger.info(f'Running server {self.__class__.__name__} on {address}{port}')

        # Finally start the server
        self.server.start()

    @property
    def address(self):
        """
        Get the address that this server is or was provided at.
        """
        return self._address

    @property
    def port(self):
        """
        Get the port that the server is or was provided on. This is 0 if a port
        was unable to be chosen.
        """
        return self._port

    @property
    def address_and_port(self) -> Tuple[str, int]:
        """
        Gets the address and port that the server is or was provided on as a tuple.

        :return: The address and port that the server is or was provided on as a tuple.
        """
        return self.address, self.port

    def setup_services(self):
        """
        Inheritors of this class should set up any services they run.
        """
        pass

    def add_service(self, service):
        """
        Add a gRPC service to this server.

        :param service: The gRPC service to add, must have the method to add the gRPC service as the attribute
        add_to_server_method.
        """
        try:
            service.add_to_server_method(service, self.server)
        except AttributeError:
            raise AttributeError("Service implementation did not have the add_to_server_method "
                                 "as an attribute, cannot automatically add to gRPC server.")

    def close(self):
        """
        Stops the server.

        Inheritors of this class should override this method with routines to stop
        services that are running.
        """
        self.server.stop(grace=False)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


# TODO: Remove this superfluous function
def get_requested_port_or_default(port: Optional[int], default: int) -> int:
    """
    Returns the port you asked for, or the default one is `port` is `None`.
    """
    if port is None:
        port = default
    return port


def is_local_address(address: str) -> bool:
    """Determine whether an address or host is local.

    This checks if the supplied ``address`` is a local IP-address/hostname.

    :param address: The IP-address or hostname whose locality is
        to be checked.

    :return: A boolean indicating if the specified address/hostname is local.


    """
    # Strip the brackets from the address if they are present
    address = address.strip('[]')

    # Attempt to parse the supplied `address`
    try:
        # If `address` is a valid IP address or range thereof
        ip_address = ipaddress.ip_address(address)
    except ValueError:
        try:
            # If the input is not a valid IP address, try to resolve it as a hostname
            ip_address = ipaddress.ip_address(socket.gethostbyname(address))
        except (socket.gaierror, ValueError):
            # If `address` can't be resolved or parsed, then assume that it is
            # non-local and let it error out somewhere else in the code.
            return False

    # Check if the address is a local one
    if ((ip_address.version == 4 and ip_address.is_loopback) or
        (ip_address.version == 6 and ip_address == ipaddress.IPv6Address('::1'))):
        return True
    else:
        return False


def is_private_address(address: str) -> bool:
    """Check if an address or host is private (i.e., used for a LAN).

    This checks if the supplied ``address`` is an IP address in one of the
    ranges reserved for private networks (10.0.0.0/8, 172.16.0.0/12, or
    192.168.0.0/16 in IPv4; fd00::/8 in IPv6).

    :param address: The IP address or hostname to be checked.

    :return: A boolean indicating if the specified address/hostname is private.
    """
    # Strip the brackets from the address if they are present
    address = address.strip('[]')

    # Attempt to parse the supplied `address`
    try:
        # If `address` is a valid IP address or range thereof
        ip_address = ipaddress.ip_address(address)
    except ValueError:
        try:
            # If the input is not a valid IP address, try to resolve it as a hostname
            ip_address = ipaddress.ip_address(socket.gethostbyname(address))
        except (socket.gaierror, ValueError):
            # If `address` can't be resolved or parsed, then assume that it is
            # non-local and let it error out somewhere else in the code.
            return False

    # Check if the address is a private one
    return ip_address.is_private


def port_in_use(port, address):
    """Check if a port is in use.

    :param port: The port to check.
    :param address: The address to check the port on.

    :return: False if the port is free, True if it is in use.
    """
    # First try IPV4 then if that errors out try IPV6
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind((address, port))
            return False
    except socket.error:
        try:
            with socket.socket(socket.AF_INET6, socket.SOCK_STREAM) as s:
                s.bind((address, port))
                return False
        except socket.error:
            return True


def _network_guard(address, port, secure, suppress_warnings):
    """Perform simple network safety guard checks.

    This function just abstracts some very verbose warning messages from the
    main body of the `GrpcServer` subroutine. These checks issue warnings when
    there is a possibility of opening up an unsecure server to external network
    traffic.
    """

    if (not secure and not is_private_address(address)) and not suppress_warnings:
        # Opening up an unsecure and unencrypted server to external network
        # traffic poses a serious security risk, and is thus discouraged.
        #
        # Admittedly this is a little heavy-handed, as allowing the server to
        # listen on all accessible IP addresses is commonly not enough to open
        # up the server to external network traffic. One must also set up some
        # form of port forwarding. However, it is hard to check locally whether
        # this has been done and thus cannot be guaranteed. So, to be on the
        # safe side this warning is always presented to the user whenever an
        # unsecured server is provided with a non-private IP address.
        #
        # This warning can be disabled by issuing suppress_warnings=True
        warnings.warn(
            textwrap.dedent(
                """
                Security Warning: Attempting to launch an insecure server on a non-private
                IP address or hostname poses a serious security risk!

                Depending on how you network is configured, your server may be open to
                unauthorized access from the internet, which may lead to data breaches or
                other security incidents. This can expose sensitive data, allow unauthorised
                control of your application, or even lead to malicious activity being conducted
                from your server.
                
                Please either:
                    1) Run your server on a local or private IP address or hostname to limit
                       accessibility to your local network only.
                    2) Implement server/client-side authentication and encryption if you need
                       to expose your server to the wider internet.
                    3) If you are certain that your server is not reachable from the internet,
                       i.e. no port forwarding has been configured, then you may suppress this
                       warning by passing `suppress_warnings=True` to the server constructor.

                If you are attempting to set up a server for remote access and are handling
                sensitive data we strongly recommend consulting a security professional to
                ensure that your server and network are configured securely.
                """).strip(), UserWarning)

    elif (not secure and not is_local_address(address)) and not suppress_warnings:
        # Setting up an unsecure server on a private network has its risks.
        # Hence, a warning is issued. However, if the user trust the network
        # then the warning can be safely ignored.
        warnings.warn(
            textwrap.dedent(
                """
                Security Warning: You're setting up an unencrypted server on a private
                network. While this is acceptable in trusted, controlled environments, it
                might expose sensitive data if untrusted clients have access to your network.
                Proceed with caution, and consider enabling encryption if you are unsure about
                the security of your network. This warning can be suppressed by passing
                `suppress_warnings=True` to the server constructor.
                """).strip(), UserWarning)

    # While we are here check if the port is already in use. No need to perform
    # a port check for port == 0 as this means "find any free port".
    if port != 0 and port_in_use(port, address):
        raise ConnectionError(
            f"The server cannot bind to the port {port} as it is already in use")
