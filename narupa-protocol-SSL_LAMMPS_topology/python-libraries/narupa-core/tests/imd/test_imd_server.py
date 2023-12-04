from typing import Generator, Tuple

import pytest
from narupa.imd.imd_client import ImdClient
from narupa.imd.imd_server import ImdServer
from narupa.imd.particle_interaction import ParticleInteraction


@pytest.fixture
def imd_server() -> Generator[ImdServer, None, None]:
    with ImdServer(address='localhost', port=0) as server:
        yield server


@pytest.fixture
def imd_server_client(imd_server) -> Generator[Tuple[ImdServer, ImdClient], None, None]:
    with ImdClient.establish_channel(address='localhost', port=imd_server.port) as client:
        yield imd_server, client


@pytest.fixture
def interaction():
    return ParticleInteraction()



