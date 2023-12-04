# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: narupa/protocol/instance/instance_service.proto
"""Generated protocol buffer code."""
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import symbol_database as _symbol_database
from google.protobuf.internal import builder as _builder
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


from google.protobuf import struct_pb2 as google_dot_protobuf_dot_struct__pb2
from google.protobuf import field_mask_pb2 as google_dot_protobuf_dot_field__mask__pb2
from narupa.protocol.instance import connect_to_trajectory_pb2 as narupa_dot_protocol_dot_instance_dot_connect__to__trajectory__pb2


DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n/narupa/protocol/instance/instance_service.proto\x12\x18narupa.protocol.instance\x1a\x1cgoogle/protobuf/struct.proto\x1a google/protobuf/field_mask.proto\x1a\x34narupa/protocol/instance/connect_to_trajectory.proto\":\n\x15LoadTrajectoryRequest\x12\x13\n\x0binstance_id\x18\x01 \x01(\t\x12\x0c\n\x04path\x18\x02 \x01(\t\"-\n\x16LoadTrajectoryResponse\x12\x13\n\x0binstance_id\x18\x01 \x01(\t\"Y\n\x16GetInstanceInfoRequest\x12\x13\n\x0binstance_id\x18\x01 \x01(\t\x12*\n\x06\x66ields\x18\x02 \x01(\x0b\x32\x1a.google.protobuf.FieldMask\"U\n\x17GetInstanceInfoResponse\x12\x13\n\x0binstance_id\x18\x01 \x01(\t\x12%\n\x04info\x18\x02 \x01(\x0b\x32\x17.google.protobuf.Struct\"+\n\x14\x43loseInstanceRequest\x12\x13\n\x0binstance_id\x18\x01 \x01(\t\",\n\x15\x43loseInstanceResponse\x12\x13\n\x0binstance_id\x18\x01 \x01(\t2\xc2\x03\n\x0fInstanceService\x12s\n\x0eLoadTrajectory\x12/.narupa.protocol.instance.LoadTrajectoryRequest\x1a\x30.narupa.protocol.instance.LoadTrajectoryResponse\x12P\n\x13\x43onnectToTrajectory\x12\x1b.ConnectToTrajectoryRequest\x1a\x1c.ConnectToTrajectoryResponse\x12v\n\x0fGetInstanceInfo\x12\x30.narupa.protocol.instance.GetInstanceInfoRequest\x1a\x31.narupa.protocol.instance.GetInstanceInfoResponse\x12p\n\rCloseInstance\x12..narupa.protocol.instance.CloseInstanceRequest\x1a/.narupa.protocol.instance.CloseInstanceResponseb\x06proto3')

_globals = globals()
_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, _globals)
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, 'narupa.protocol.instance.instance_service_pb2', _globals)
if _descriptor._USE_C_DESCRIPTORS == False:
  DESCRIPTOR._options = None
  _globals['_LOADTRAJECTORYREQUEST']._serialized_start=195
  _globals['_LOADTRAJECTORYREQUEST']._serialized_end=253
  _globals['_LOADTRAJECTORYRESPONSE']._serialized_start=255
  _globals['_LOADTRAJECTORYRESPONSE']._serialized_end=300
  _globals['_GETINSTANCEINFOREQUEST']._serialized_start=302
  _globals['_GETINSTANCEINFOREQUEST']._serialized_end=391
  _globals['_GETINSTANCEINFORESPONSE']._serialized_start=393
  _globals['_GETINSTANCEINFORESPONSE']._serialized_end=478
  _globals['_CLOSEINSTANCEREQUEST']._serialized_start=480
  _globals['_CLOSEINSTANCEREQUEST']._serialized_end=523
  _globals['_CLOSEINSTANCERESPONSE']._serialized_start=525
  _globals['_CLOSEINSTANCERESPONSE']._serialized_end=569
  _globals['_INSTANCESERVICE']._serialized_start=572
  _globals['_INSTANCESERVICE']._serialized_end=1022
# @@protoc_insertion_point(module_scope)