# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Optional
from narupa.core import NarupaServer, get_requested_port_or_default, DEFAULT_SERVE_ADDRESS, GrpcCredentials
from narupa.protocol.trajectory import add_TrajectoryServiceServicer_to_server
from .frame_data import FrameData
from .frame_publisher import FramePublisher

DEFAULT_PORT = 54321

PLAY_COMMAND_KEY = "playback/play"
RESET_COMMAND_KEY = "playback/reset"
STEP_COMMAND_KEY = "playback/step"
PAUSE_COMMAND_KEY = "playback/pause"
GET_DYNAMICS_INTERVAL_COMMAND_KEY = 'trajectory/get-dynamics-interval'
SET_DYNAMICS_INTERVAL_COMMAND_KEY = 'trajectory/set-dynamics-interval'


class FrameServer(NarupaServer):
    _trajectory_service: FramePublisher

    def __init__(
            self, *, address: Optional[str] = None, port: Optional[int] = None,
            credentials: Optional[GrpcCredentials] = None,):
        if address is None:
            address = DEFAULT_SERVE_ADDRESS
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        super().__init__(address=address, port=port, credentials=credentials)
        self._frame_count = 0

    def setup_services(self):
        super().setup_services()
        self._trajectory_service = FramePublisher()
        self._trajectory_service.add_to_server_method(self._trajectory_service, self.server)

    def send_frame(self, frame_index: int, frame_data: FrameData):
        self._trajectory_service.send_frame(frame_index, frame_data.raw)
        self._frame_count += 1

    def close(self):
        super().close()
        self._trajectory_service.close()

    @property
    def frame_count(self):
        """
        Counts how many times send_frame has been called on this publisher.
        """
        return self._frame_count
