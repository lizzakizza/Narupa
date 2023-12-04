# Generated by the gRPC Python protocol compiler plugin. DO NOT EDIT!
"""Client and server classes corresponding to protobuf-defined services."""
import grpc

from narupa.protocol.command import command_service_pb2 as narupa_dot_protocol_dot_command_dot_command__service__pb2


class CommandStub(object):
    """Missing associated documentation comment in .proto file."""

    def __init__(self, channel):
        """Constructor.

        Args:
            channel: A grpc.Channel.
        """
        self.GetCommands = channel.unary_unary(
                '/narupa.protocol.command.Command/GetCommands',
                request_serializer=narupa_dot_protocol_dot_command_dot_command__service__pb2.GetCommandsRequest.SerializeToString,
                response_deserializer=narupa_dot_protocol_dot_command_dot_command__service__pb2.GetCommandsReply.FromString,
                )
        self.RunCommand = channel.unary_unary(
                '/narupa.protocol.command.Command/RunCommand',
                request_serializer=narupa_dot_protocol_dot_command_dot_command__service__pb2.CommandMessage.SerializeToString,
                response_deserializer=narupa_dot_protocol_dot_command_dot_command__service__pb2.CommandReply.FromString,
                )


class CommandServicer(object):
    """Missing associated documentation comment in .proto file."""

    def GetCommands(self, request, context):
        """Get a list of all the commands available on this service 
        """
        context.set_code(grpc.StatusCode.UNIMPLEMENTED)
        context.set_details('Method not implemented!')
        raise NotImplementedError('Method not implemented!')

    def RunCommand(self, request, context):
        """Runs a command on the service 
        """
        context.set_code(grpc.StatusCode.UNIMPLEMENTED)
        context.set_details('Method not implemented!')
        raise NotImplementedError('Method not implemented!')


def add_CommandServicer_to_server(servicer, server):
    rpc_method_handlers = {
            'GetCommands': grpc.unary_unary_rpc_method_handler(
                    servicer.GetCommands,
                    request_deserializer=narupa_dot_protocol_dot_command_dot_command__service__pb2.GetCommandsRequest.FromString,
                    response_serializer=narupa_dot_protocol_dot_command_dot_command__service__pb2.GetCommandsReply.SerializeToString,
            ),
            'RunCommand': grpc.unary_unary_rpc_method_handler(
                    servicer.RunCommand,
                    request_deserializer=narupa_dot_protocol_dot_command_dot_command__service__pb2.CommandMessage.FromString,
                    response_serializer=narupa_dot_protocol_dot_command_dot_command__service__pb2.CommandReply.SerializeToString,
            ),
    }
    generic_handler = grpc.method_handlers_generic_handler(
            'narupa.protocol.command.Command', rpc_method_handlers)
    server.add_generic_rpc_handlers((generic_handler,))


 # This class is part of an EXPERIMENTAL API.
class Command(object):
    """Missing associated documentation comment in .proto file."""

    @staticmethod
    def GetCommands(request,
            target,
            options=(),
            channel_credentials=None,
            call_credentials=None,
            insecure=False,
            compression=None,
            wait_for_ready=None,
            timeout=None,
            metadata=None):
        return grpc.experimental.unary_unary(request, target, '/narupa.protocol.command.Command/GetCommands',
            narupa_dot_protocol_dot_command_dot_command__service__pb2.GetCommandsRequest.SerializeToString,
            narupa_dot_protocol_dot_command_dot_command__service__pb2.GetCommandsReply.FromString,
            options, channel_credentials,
            insecure, call_credentials, compression, wait_for_ready, timeout, metadata)

    @staticmethod
    def RunCommand(request,
            target,
            options=(),
            channel_credentials=None,
            call_credentials=None,
            insecure=False,
            compression=None,
            wait_for_ready=None,
            timeout=None,
            metadata=None):
        return grpc.experimental.unary_unary(request, target, '/narupa.protocol.command.Command/RunCommand',
            narupa_dot_protocol_dot_command_dot_command__service__pb2.CommandMessage.SerializeToString,
            narupa_dot_protocol_dot_command_dot_command__service__pb2.CommandReply.FromString,
            options, channel_credentials,
            insecure, call_credentials, compression, wait_for_ready, timeout, metadata)
